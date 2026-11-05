"""Direct negative-synergy learning, calibration and ranking evaluation."""
from dataclasses import dataclass
from pathlib import Path
import json
import time
import warnings
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from scipy.special import expit
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss, log_loss
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from .training import FeatureNormalizer
from .models import invariant_pair, normalized_graph
from .residual_models import GatedGeneEncoder, ResidualModelConfig


def load_inputs(root):
    labels = pd.read_csv(root/'data/sl/ht29_labels.csv')
    labels = labels[labels.used_for_model].reset_index(drop=True)
    genes = pd.read_csv(root/'data/processed/model_genes.csv').gene.tolist()
    lookup = {g:i for i,g in enumerate(genes)}
    pairs = labels[['gene_a','gene_b']].apply(lambda c:c.map(lookup)).to_numpy(np.int64)
    features = np.load(root/'data/evidence/integrated_features.npz')['features']
    graphs = dict(np.load(root/'data/evidence/graphs.npz'))
    ann = pd.read_csv(root/'data/pairs_sources/HT29_annotations.csv')
    mapping = pd.read_csv(root/'data/annotations/symbol_mapping.csv').set_index('source_symbol').gene.to_dict()
    anchor = set(ann.loc[ann.MyNote.eq('AnchorSingletons'),'Gene2'].map(lambda g:mapping.get(g,g)))
    return labels, genes, pairs, features, graphs, anchor


def partition(table, protocol, config, anchors):
    seed = config['partition_seed']; rng = np.random.default_rng(seed)
    fractions=config.get('fractions',{'validation':.15,'calibration':.1,'test':.15})
    fv,fc,ft=(float(fractions[k]) for k in ['validation','calibration','test'])
    if min(fv,fc,ft)<=0 or fv+fc+ft>=1:
        raise ValueError('Validation, calibration and test fractions must be positive and sum to less than one.')
    ids = np.arange(len(table))
    if protocol=='pair_holdout':
        # Stratification is fixed once, shared across every architecture/seed.
        train, held = train_test_split(ids,test_size=fv+fc+ft,random_state=seed,stratify=table.synergy_label)
        val, rest = train_test_split(held,test_size=(fc+ft)/(fv+fc+ft),random_state=seed,stratify=table.synergy_label.iloc[held])
        cal, test = train_test_split(rest,test_size=ft/(fc+ft),random_state=seed,stratify=table.synergy_label.iloc[rest])
        return dict(train=train,val=val,calibration=cal,test=test,discarded=np.array([],int))
    universe = np.array(sorted(set(table.gene_a)|set(table.gene_b)))
    if protocol=='anchor_holdout':
        universe = np.array(sorted(set(universe)&anchors))
    universe = rng.permutation(universe)
    n=len(universe); nt=round(n*ft); nc=round(n*fc); nv=round(n*fv)
    groups={'test':set(universe[:nt]),'calibration':set(universe[nt:nt+nc]),'val':set(universe[nt+nc:nt+nc+nv])}
    membership={name:table.gene_a.isin(gs)|table.gene_b.isin(gs) for name,gs in groups.items()}
    nheld=sum(membership.values()); out={'train':ids[nheld.eq(0)]}
    for name,mask in membership.items(): out[name]=ids[mask & nheld.eq(1)]
    out['discarded']=ids[nheld.gt(1)]
    return out


class SynergyNet(nn.Module):
    def __init__(self, n_features, variant, config):
        super().__init__(); self.variant=variant
        self.graph_names={'string':('string',),'multigraph':('reactome','string','bioplex')}.get(variant,())
        if variant=='linear':
            self.encoder=nn.Identity(); self.head=nn.Linear(3*n_features,3)
        else:
            h=config['hidden_dim']
            self.encoder=GatedGeneEncoder(n_features,ResidualModelConfig(hidden_dim=h,dropout=config['dropout'],graph_names=self.graph_names))
            self.head=nn.Sequential(nn.Linear(3*h,2*h),nn.GELU(),nn.Dropout(config['dropout']),nn.Linear(2*h,h),nn.GELU(),nn.Dropout(config['dropout']),nn.Linear(h,3))
    def encode(self,x,graphs):
        return self.encoder(x) if self.variant=='linear' else self.encoder(x,graphs)
    def forward(self,h,pairs):
        return self.head(invariant_pair(h[pairs[:,0]],h[pairs[:,1]]))


@dataclass
class Fit:
    model: object
    normalizer: object
    center: float
    scale: float
    variant: str
    device: str
    history: list
    elapsed: float
    best_epoch: int
    config: dict


def _tensors(features,pairs,normalizer,graphs,names,device):
    x=torch.as_tensor(normalizer.transform(features),device=device)
    p=torch.as_tensor(pairs,device=device)
    gs={g:normalized_graph(graphs[g],len(features),device) for g in names}
    return x,p,gs


def fit(features,pairs,labels,graphs,train,val,variant,seed,config,device):
    torch.set_num_threads(4); torch.manual_seed(seed); np.random.seed(seed)
    normalizer=FeatureNormalizer.fit(features,np.unique(pairs[train]))
    delta=labels.dlfc.to_numpy(np.float32)
    center=float(np.mean(delta[train])); scale=max(float(np.std(delta[train])),.05)
    start=time.monotonic()
    if variant=='xgb':
        import xgboost as xgb
        x=normalizer.transform(features); a,b=x[pairs[:,0]],x[pairs[:,1]]
        z=np.concatenate([a+b,np.abs(a-b),a*b],axis=1)
        ds=xgb.QuantileDMatrix(z[train],label=(delta[train]<-1).astype(float),max_bin=64)
        dv=xgb.QuantileDMatrix(z[val],label=(delta[val]<-1).astype(float),ref=ds,max_bin=64)
        both_classes=len(np.unique(delta[val]<-1))==2
        params=dict(device=device,tree_method='hist',max_depth=3,max_bin=64,eta=.04,subsample=.8,colsample_bytree=.8,reg_lambda=10,nthread=4,seed=seed,objective='binary:logistic',eval_metric='aucpr' if both_classes else 'logloss')
        h={}; classifier=xgb.train(params,ds,500,[(dv,'validation')],early_stopping_rounds=25,maximize=both_classes,evals_result=h,verbose_eval=False)
        ds.set_label(delta[train]);dv.set_label(delta[val]);params.update(objective='reg:squarederror',eval_metric='rmse')
        reg=xgb.train(params,ds,500,[(dv,'validation')],early_stopping_rounds=25,verbose_eval=False)
        return Fit((classifier,reg),normalizer,center,scale,variant,device,[],time.monotonic()-start,classifier.best_iteration+1,config)
    model=SynergyNet(features.shape[1],variant,config).to(device)
    x,p,gs=_tensors(features,pairs,normalizer,graphs,model.graph_names,device)
    y=torch.as_tensor(np.column_stack([delta<-1,delta<-.75,(delta-center)/scale]),dtype=torch.float32,device=device)
    ti=torch.as_tensor(train,device=device);vi=torch.as_tensor(val,device=device)
    positives=y[ti,:2].sum(0); pw=((len(train)-positives)/(positives+1)).sqrt().clamp(1,20)
    opt=torch.optim.AdamW(model.parameters(),lr=config['learning_rate'],weight_decay=config['weight_decay'])
    best=-np.inf;best_state=None;stale=0;history=[];best_epoch=0
    for epoch in range(1,config['max_epochs']+1):
        model.train();total=0
        for batch in ti[torch.randperm(len(ti),device=device)].split(config['batch_size']):
            opt.zero_grad(set_to_none=True);out=model(model.encode(x,gs),p[batch]);target=y[batch]
            regression=F.smooth_l1_loss(out[:,2],target[:,2])
            if variant=='regression_only': loss=regression
            else:
                classification=F.binary_cross_entropy_with_logits(out[:,:2],target[:,:2],pos_weight=pw)
                pos=out[target[:,0].eq(1),0];neg=out[target[:,0].eq(0),0]
                ranking=F.softplus(neg[:512,None]-pos[None,:]).mean() if len(pos) and len(neg) else out.sum()*0
                loss=classification+config['regression_weight']*regression+config['ranking_weight']*ranking
            loss.backward();nn.utils.clip_grad_norm_(model.parameters(),5);opt.step();total+=float(loss.detach())*len(batch)
        model.eval()
        with torch.no_grad(): vout=model(model.encode(x,gs),p[vi]).cpu().numpy()
        truth=delta[val]<-1;score=-vout[:,2] if variant=='regression_only' else vout[:,0]
        rmse=float(np.sqrt(np.mean((vout[:,2]*scale+center-delta[val])**2)))
        ap=float(average_precision_score(truth,score)) if truth.any() else np.nan
        criterion=-rmse if variant=='regression_only' else (ap if np.isfinite(ap) else -float(np.mean(np.logaddexp(0,score)-truth*score)))
        history.append(dict(epoch=epoch,train_loss=total/len(train),validation_ap=ap,validation_rmse=rmse))
        if criterion>best+1e-7:
            best=criterion;best_epoch=epoch;stale=0;best_state={k:v.detach().clone() for k,v in model.state_dict().items()}
        else:stale+=1
        if stale>=config['patience']:break
    model.load_state_dict(best_state);model.eval()
    return Fit(model,normalizer,center,scale,variant,device,history,time.monotonic()-start,best_epoch,config)


def predict(result,features,pairs,graphs,passes=1,seed=0):
    if result.variant=='xgb':
        import xgboost as xgb
        x=result.normalizer.transform(features);a,b=x[pairs[:,0]],x[pairs[:,1]]
        dm=xgb.DMatrix(np.concatenate([a+b,np.abs(a-b),a*b],axis=1))
        clf,reg=result.model
        score=clf.predict(dm,output_margin=True,iteration_range=(0,clf.best_iteration+1))
        dlfc=reg.predict(dm,iteration_range=(0,reg.best_iteration+1))
        return dict(score=score,dlfc=dlfc,uncertainty=np.zeros(len(pairs)))
    torch.manual_seed(seed);model=result.model
    x,p,gs=_tensors(features,pairs,result.normalizer,graphs,model.graph_names,result.device)
    model.train(passes>1);scores=[];deltas=[]
    with torch.no_grad():
        for _ in range(passes):
            h=model.encode(x,gs);out=torch.cat([model(h,batch) for batch in p.split(4096)]).cpu().numpy()
            deltas.append(out[:,2]*result.scale+result.center)
            scores.append(-out[:,2] if result.variant=='regression_only' else out[:,0])
    model.eval();s=np.stack(scores)
    return dict(score=s.mean(0),dlfc=np.mean(deltas,axis=0),uncertainty=expit(s).std(0))


def calibrate(scores,truth):
    truth=np.asarray(truth,int);scores=np.asarray(scores)
    if len(np.unique(truth))<2:
        return {'slope':0.,'intercept':float(np.log((truth.sum()+.5)/(len(truth)-truth.sum()+.5))),'status':'constant_prevalence_no_two_classes','n':len(truth),'positives':int(truth.sum())}
    lr=LogisticRegression(C=1.,solver='lbfgs').fit(scores.reshape(-1,1),truth)
    # Reversed calibration indicates a failed score, never silently invert ranking.
    slope=max(0.,float(lr.coef_[0,0]));intercept=float(lr.intercept_[0])
    if slope==0:intercept=float(np.log((truth.sum()+.5)/(len(truth)-truth.sum()+.5)))
    return {'slope':slope,'intercept':intercept,'status':'platt_positive_slope' if slope>0 else 'constant_prevalence_reversed_fit','n':len(truth),'positives':int(truth.sum())}


def probability(cal,scores):return expit(cal['slope']*np.asarray(scores)+cal['intercept'])


def metrics(table,score,delta=None,prob=None):
    y=table.synergy_label.to_numpy(int);n=len(y);pos=int(y.sum());order=np.argsort(-np.asarray(score),kind='stable')
    out={'n':n,'positives':pos,'prevalence':float(y.mean()),'auprc':float(average_precision_score(y,score)) if pos else np.nan,
         'auroc':float(roc_auc_score(y,score)) if 0<pos<n else np.nan}
    for k in [20,50,100,500]:
        kk=min(k,n);hits=int(y[order[:kk]].sum());out.update({f'precision_at_{k}':hits/kk,f'recall_at_{k}':hits/pos if pos else np.nan,f'enrichment_at_{k}':hits/kk/y.mean() if pos else np.nan,f'hits_at_{k}':hits})
    if delta is not None:
        out['dlfc_rmse']=float(np.sqrt(np.mean((table.dlfc.to_numpy()-delta)**2)))
        out['dlfc_spearman']=float(spearmanr(table.dlfc,delta).statistic)
    if prob is not None:
        out['brier']=float(brier_score_loss(y,prob));out['logloss']=float(log_loss(y,prob,labels=[0,1]))
    for t in [.75,1,1.25]:
        target=table.dlfc.lt(-t).to_numpy(int)
        out[f'auprc_dlfc_lt_{t}']=float(average_precision_score(target,score)) if target.sum() else np.nan
    return out


def save(result,directory):
    directory.mkdir(parents=True,exist_ok=True)
    torch.save(result,directory/'fit.pt')
    pd.DataFrame(result.history).to_csv(directory/'history.csv',index=False)


def load(directory,device):
    result=torch.load(directory/'fit.pt',map_location=device,weights_only=False);result.device=device
    if result.variant=='xgb':
        for model in result.model:model.set_param({'device':device})
    return result
