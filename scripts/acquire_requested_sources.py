"""Download fixed public inputs; preserve URLs, bytes and SHA256 per file.

Uses the Python standard library so acquisition can run on either host.
The numerical processing step is prepare_integrated_inputs.py.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request
import csv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/evidence/requested_sources'
FILES = {
    'DepMap_Model.csv': 'https://ndownloader.figshare.com/files/51065297',
    'DepMap_CRISPRGeneEffect.csv': 'https://ndownloader.figshare.com/files/51064667',
    'TCGA_COAD_expression.tsv.gz': 'https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.gz',
    'TCGA_COAD_metadata.json': 'https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.json',
    'GSE132465_cell_annotation.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz',
    'GSE132465_counts.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz',
    'DepMap_release.json': 'https://api.figshare.com/v2/articles/27993248',
    'HCT116_flux_Control_clean.p': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Results/pTVA_Control_clean.p',
    'HCT116_flux_OxR_clean.p': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Results/pTVA_OxR_clean.p',
    'HCT116_metabolic_model.mat': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Models/HCT_GEM.mat',
    'HCT116_flux_analysis.ipynb': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Code/flux_analysis.ipynb',
    'HCT116_exchange_constraints.ipynb': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Code/exometabo_constrain.ipynb',
    'HCT116_growth_curves.csv': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Data/growthCurves.csv',
    'HCT116_metabolite_time_series.csv': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Data/metabolites_nM_timeSeriesAnnotated.csv',
    'hg19_refGene.txt.gz': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz',
    'GSM4115600_H3K4me3_HCT116_peaks.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4115nnn/GSM4115600/suppl/GSM4115600_H3K4me3_HCT116_peaks.txt.gz',
    'GSM4115601_H3K4me3_HCT116OxR_peaks.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4115nnn/GSM4115601/suppl/GSM4115601_H3K4me3_HCT116OxR_peaks.txt.gz',
    'GSM4115602_H3K27ac_HCT116_peaks.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4115nnn/GSM4115602/suppl/GSM4115602_H3K27ac_HCT116_peaks.txt.gz',
    'GSM4115603_H3K27ac_HCT116OxR_peaks.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4115nnn/GSM4115603/suppl/GSM4115603_H3K27ac_HCT116OxR_peaks.txt.gz',
    'GSM4115604_H3K27me3_HCT116_peaks.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4115nnn/GSM4115604/suppl/GSM4115604_H3K27me3_HCT116_peaks.txt.gz',
    'GSM4115605_H3K27me3_HCT116OxR_peaks.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4115nnn/GSM4115605/suppl/GSM4115605_H3K27me3_HCT116OxR_peaks.txt.gz',
    'GSM4115600_metadata.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4115600&targ=self&form=text&view=full',
    'HCT116_general_functions.ipynb': 'https://raw.githubusercontent.com/HAHerrmann/Hct116_DrugRes/3ba911ee3b931e0909add599beba48bdca9f69a9/Code/general_functions.ipynb',
}


def segmented(name, url, size):
    """Bounded range reads; resume complete pieces and verify server ranges."""
    path = OUT / name
    pieces = OUT / (name + '.segments')
    pieces.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + '.part')
    prefix = partial.stat().st_size if partial.exists() else 0
    chunk = 4 * 2**20
    def read_piece(start):
        end = min(size, start + chunk)
        piece = pieces / str(start)
        if piece.exists() and piece.stat().st_size == end - start:
            return piece
        if end <= prefix:
            with partial.open('rb') as stream:
                stream.seek(start)
                piece.write_bytes(stream.read(end - start))
            return piece
        for attempt in range(4):
            try:
                request = urllib.request.Request(url, headers={'Range': f'bytes={start}-{end-1}'})
                with urllib.request.urlopen(request, timeout=90) as response:
                    expected = f'bytes {start}-{end-1}/{size}'
                    if response.status != 206 or response.headers.get('Content-Range') != expected:
                        raise ValueError('Server did not honor the requested byte range')
                    content = response.read()
                    if len(content) != end-start:
                        raise ValueError('Incomplete piece')
                    piece.write_bytes(content)
                return piece
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(attempt + 1)
    starts = list(range(0, size, chunk))
    with ThreadPoolExecutor(max_workers=10) as executor:
        for number, future in enumerate(as_completed([executor.submit(read_piece, start) for start in starts]), start=1):
            future.result()
            if number % 10 == 0:
                print(json.dumps({'file': name, 'completed_pieces': number, 'total_pieces': len(starts)}), flush=True)
    with path.open('wb') as stream:
        for start in starts:
            stream.write((pieces / str(start)).read_bytes())
    import shutil
    shutil.rmtree(pieces)
    if partial.exists():
        partial.unlink()


def acquire(name, url, data=None):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    sidecar = path.with_name(path.name + '.provenance.json')
    if path.exists() and sidecar.exists():
        return json.loads(sidecar.read_text())
    known_sizes = {'DepMap_CRISPRGeneEffect.csv': 428678699, 'GSE132465_counts.txt.gz': 128941455}
    if name in known_sizes and data is None:
        segmented(name, url, known_sizes[name])
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(2**20), b''):
                digest.update(block)
        record = {'file': name, 'url': url, 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size,
                  'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'status': 'downloaded'}
        sidecar.write_text(json.dumps(record, indent=2) + '\n')
        print(json.dumps(record), flush=True)
        return record
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, data=data, headers={'User-Agent': 'MetaGene academic research contact via public APIs'})
            digest = hashlib.sha256()
            with urllib.request.urlopen(request, timeout=60) as response, path.with_suffix(path.suffix + '.part').open('wb') as stream:
                while True:
                    block = response.read(2**20)
                    if not block:
                        break
                    digest.update(block)
                    stream.write(block)
            path.with_suffix(path.suffix + '.part').replace(path)
            record = {'file': name, 'url': url, 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size,
                      'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'status': 'downloaded'}
            sidecar.write_text(json.dumps(record, indent=2) + '\n')
            print(json.dumps(record), flush=True)
            return record
        except Exception as exc:
            print(json.dumps({'file': name, 'attempt': attempt + 1, 'error': str(exc)}), flush=True)
            time.sleep(2 * (attempt + 1))
    return {'file': name, 'url': url, 'status': 'unavailable_after_three_attempts'}


def string_inputs():
    genes = list(csv.DictReader((ROOT / 'data/processed/model_genes.csv').open()))
    identifiers = '\r'.join(row['gene'] for row in genes)
    params = {'identifiers': identifiers, 'species': 9606, 'echo_query': 1, 'caller_identity': 'MetaGene'}
    base = 'https://version-12-5.string-db.org/api/tsv/'
    mapping = acquire('STRING_mapping.tsv', base + 'get_string_ids', urllib.parse.urlencode(params).encode())
    if mapping['status'] != 'downloaded':
        return [mapping]
    time.sleep(1)
    rows = list(csv.DictReader((OUT / 'STRING_mapping.tsv').open(), delimiter='\t'))
    params = {'identifiers': '\r'.join(row['stringId'] for row in rows), 'species': 9606,
              'required_score': 700, 'network_type': 'physical', 'add_nodes': 0, 'caller_identity': 'MetaGene'}
    network = acquire('STRING_physical_network.tsv', base + 'network', urllib.parse.urlencode(params).encode())
    return [mapping, network]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', nargs='*', choices=list(FILES) + ['STRING'])
    args = parser.parse_args()
    chosen = args.only or list(FILES) + ['STRING']
    records = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        jobs = [executor.submit(acquire, name, url) for name, url in FILES.items() if name in chosen]
        for job in as_completed(jobs):
            records.append(job.result())
    if 'STRING' in chosen:
        records.extend(string_inputs())
    all_records = {json.loads(path.read_text())['file']: json.loads(path.read_text())
                   for path in OUT.glob('*.provenance.json')}
    all_records.update({record['file']: record for record in records})
    (OUT / 'acquisition_manifest.json').write_text(json.dumps(
        [all_records[name] for name in sorted(all_records)], indent=2) + '\n')


if __name__ == '__main__':
    main()
