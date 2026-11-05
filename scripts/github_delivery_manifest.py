"""Describe the complete GitHub project without creating a delivery archive."""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = 'https://github.com/Tairan-Terrian/MetaGene'


def main():
    manifest = ROOT / 'results/delivery_manifest.json'
    listed = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=ROOT)
    files = sorted({name.decode('utf-8') for name in listed.split(b'\0') if name})
    records = []
    for relative in files:
        path = ROOT / relative
        if path == manifest:
            continue
        if not path.is_file():
            raise FileNotFoundError(path)
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(2**20), b''):
                digest.update(block)
        records.append({'path': relative, 'bytes': path.stat().st_size,
                        'sha256': digest.hexdigest()})
    payload = {
        'delivery_type': 'Complete GitHub project with Git LFS',
        'repository': REPOSITORY,
        'version': json.loads((ROOT / 'configs/experiment_sl.json').read_text())['version'],
        'file_count': len(records),
        'total_uncompressed_bytes': sum(row['bytes'] for row in records),
        'scope': 'Implementation code, configuration, original data snapshots, processed inputs, model weights, per-run outputs, reports, figures, candidate tables and provenance. Environments, caches, new local reproduction runs, local-only bookkeeping and the retired delivery archive are excluded by .gitignore.',
        'hash_basis': 'SHA256 of original working-tree bytes, including expanded Git LFS content. The manifest itself is excluded from its own hashes. The enclosing Git commit identifies this delivery snapshot.',
        'files': records,
    }
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({key: payload[key] for key in ['repository', 'file_count', 'total_uncompressed_bytes']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
