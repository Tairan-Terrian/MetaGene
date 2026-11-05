"""Render the current synergistic-lethality report from completed measured results."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/sl'
COLORS={'linear':'#767676','xgb':'#B64342','graphless':'#AADCA9','string':'#3775BA','multigraph':'#0F4D92','regression_only':'#E9A6A1'}
NAMES={'linear':'Linear multitask','xgb':'XGBoost','graphless':'No graph','string':'STRING GNN','multigraph':'Multigraph GNN','regression_only':'Regression only'}


def md_table(frame,columns=None):
    if columns:frame=frame[columns]
    frame=frame.copy()
    for c in frame:
        if pd.api.types.is_float_dtype(frame[c]):frame[c]=frame[c].map(lambda x:'不可评价' if pd.isna(x) else f'{x:.4f}')
    return '| '+' | '.join(frame.columns)+' |\n|'+'|'.join(['---']*len(frame.columns))+'|\n'+'\n'.join('| '+' | '.join(str(x).replace('|',' / ') for x in row)+' |' for row in frame.itertuples(index=False,name=None))


def figures(selected):
    figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':1.5,'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42})
    def save(fig,name):
        fig.tight_layout(pad=1.5)
        for ext in ['png','pdf','svg']:fig.savefig(figdir/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    labels=pd.read_csv(ROOT/'data/sl/ht29_labels.csv');labels=labels[labels.used_for_model]
    fig,axes=plt.subplots(1,2,figsize=(13,4.8))
    axes[0].hist(labels.dlfc,bins=70,color='#3775BA',edgecolor='white',linewidth=.3)
    axes[0].axvline(-1,color='#B64342',ls='--',label='Strong synergy threshold: -1');axes[0].set(xlabel='Control-centered dLFC',ylabel='Measured pairs',title=f'A  Label distribution (n={len(labels):,})');axes[0].legend(fontsize=10)
    x=labels.sort_values('dlfc').head(35)
    axes[1].errorbar(np.arange(len(x)),x.dlfc,yerr=x.replicate_se,fmt='o',ms=4,color='#0F4D92',elinewidth=1,capsize=2)
    axes[1].axhline(-1,color='#B64342',ls='--');axes[1].set(xlabel='Measured pairs sorted by dLFC (first 35)',ylabel='dLFC ± replicate SE',title='B  Strong-signal tail (SE is not full assay error)')
    save(fig,'01_labels')
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for model in ['linear','xgb','regression_only',selected]:
        t=pd.read_csv(OUT/f'ensemble_pair_holdout_{model}_test.csv');p,r,_=precision_recall_curve(t.synergy_label,t.model_score)
        axes[0].step(r,p,where='post',label=NAMES[model],color=COLORS[model],lw=2)
        order=t.sort_values('model_score',ascending=False);hits=order.synergy_label.cumsum().to_numpy()
        axes[1].plot(np.arange(1,len(t)+1),hits,label=NAMES[model],color=COLORS[model],lw=2)
    axes[0].axhline(t.synergy_label.mean(),color='#999999',ls=':',label='Prevalence');axes[0].set(xlim=(0,1),ylim=(0,1.02),xlabel='Recall',ylabel='Precision',title='A  Fixed test set: 5 positives / 2,229 pairs');axes[0].legend(fontsize=9)
    axes[1].plot([0,len(t)],[0,t.synergy_label.sum()],color='#999999',ls=':',label='Random expectation');axes[1].axvline(500,color='#CFCECE',ls='--');axes[1].set(xlabel='Number selected',ylabel='Strong-synergy hits recovered',title='B  Cumulative recovery');axes[1].legend(fontsize=9)
    save(fig,'02_ranking')
    metrics=pd.read_csv(OUT/'ensemble_metrics.csv');fig,axes=plt.subplots(1,3,figsize=(16,4.8))
    for ax,protocol,title in zip(axes,['pair_holdout','gene_holdout','anchor_holdout'],['Pair holdout','Unseen gene','Unseen anchor']):
        f=metrics[metrics.protocol.eq(protocol)].set_index('model').loc[list(NAMES)]
        ax.barh(np.arange(len(f)),f.auprc,color=[COLORS[m] for m in f.index],edgecolor='#272727',linewidth=.5)
        ax.set_yticks(np.arange(len(f)),[NAMES[m] for m in f.index],fontsize=10);ax.invert_yaxis();ax.axvline(f.prevalence.iloc[0],color='#B64342',ls=':')
        ax.set(xlabel='Average precision',title=f'{title}\n{int(f.positives.iloc[0])} positives / {int(f.n.iloc[0]):,} pairs')
    save(fig,'03_benchmarks')
    ac=pd.read_csv(OUT/'active_learning_summary.csv');fig,axes=plt.subplots(1,2,figsize=(13,4.8));palette={'random':'#767676','probability':'#0F4D92','uncertainty':'#B64342','upper_confidence':'#42949E'}
    for strategy,f in ac.groupby('strategy',sort=False):
        axes[0].plot(f.n_measured,f.strong_hits_mean,'o-',ms=3,lw=2,color=palette[strategy],label=strategy)
        axes[0].fill_between(f.n_measured,f.strong_hits_mean-f.strong_hits_sd,f.strong_hits_mean+f.strong_hits_sd,alpha=.13,color=palette[strategy])
        axes[1].plot(f.n_measured,f.tail_hits_mean,'o-',ms=3,lw=2,color=palette[strategy],label=strategy)
    axes[0].set(xlabel='Cumulative measured pairs',ylabel='Hits with dLFC < -1',title='A  Strong-synergy discovery (mean ± seed SD)');axes[0].legend(fontsize=9)
    axes[1].set(xlabel='Cumulative measured pairs',ylabel='Hits with dLFC < -0.75',title='B  Auxiliary negative-tail sensitivity');axes[1].legend(fontsize=9)
    save(fig,'04_active_learning')
    top=pd.read_csv(OUT/'candidates/top500.csv');fig,axes=plt.subplots(1,2,figsize=(13,4.8))
    for tier,color,label in [(0,'#0F4D92','Singles pass initial cutoff'),(1,'#42949E','Missing single-gene evidence'),(2,'#B64342','Strong single-gene depletion')]:
        f=top[top.single_effect_tier.eq(tier)];axes[0].scatter(f.pred_dlfc,f.synergy_probability,s=15,alpha=.65,color=color,label=f'{label} (n={len(f)})')
    axes[0].axvline(-1,color='#767676',ls='--');axes[0].set(xlabel='Predicted dLFC',ylabel='Calibrated HT29 event score',title='A  Top500 prediction and eligibility');axes[0].legend(fontsize=8)
    hist=top.chromatin_gene.value_counts().sort_values();axes[1].barh(hist.index,hist.values,color='#3775BA');axes[1].set(xlabel='Pairs in Top500',title='B  Histone-enzyme composition')
    save(fig,'05_candidates')


def main():
    import argparse
    from write_sl_documents import main as write_documents
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--documents-only',action='store_true',help='Update reading documents using compact results; keep existing figures.')
    args=parser.parse_args()
    if not args.documents_only:
        selected=json.loads((OUT/'model_selection.json').read_text())['selected_gnn']
        figures(selected)
    write_documents()


if __name__=='__main__':main()
