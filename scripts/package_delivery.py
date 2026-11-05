"""Package the complete portable project, excluding caches and local environments."""
import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    if (ROOT / 'results/local_storage_manifest.json').exists():
        raise SystemExit('This is the lightweight local review copy. Package the complete project on the compute server instead.')
    output = ROOT / 'MetaGene交接包.tar.gz'
    manifest = ROOT / 'results/delivery_manifest.json'
    excluded_parts = {'.venv', '.Rlib', '.toolchain', '__pycache__', '.git', '.DS_Store'}
    files = sorted(p for p in ROOT.rglob('*') if p.is_file()
                   and not (set(p.relative_to(ROOT).parts) & excluded_parts)
                   and p not in {output, manifest}
                   and not p.name.endswith(('.pyc', '.part')))
    records = []
    for path in files:
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(2**20), b''):
                digest.update(block)
        records.append({'path': path.relative_to(ROOT).as_posix(), 'bytes': path.stat().st_size,
                        'sha256': digest.hexdigest()})
    manifest.write_text(json.dumps({'package': output.name, 'file_count': len(records),
        'total_uncompressed_bytes': sum(r['bytes'] for r in records),
        'excluded': sorted(excluded_parts),
        'scope': 'Current code, configuration, source snapshots, processed inputs, all numerical results, fitted weights and reports. The manifest itself is included but not recursively hashed.',
        'files': records}, ensure_ascii=False, indent=2) + '\n')
    with tarfile.open(output, 'w:gz', compresslevel=3) as archive:
        for path in files + [manifest]:
            archive.add(path, arcname=str(Path('MetaGene') / path.relative_to(ROOT)), recursive=False)
    expected = {str(Path('MetaGene') / r['path']): r for r in records}
    manifest_name = str(Path('MetaGene') / manifest.relative_to(ROOT))
    expected[manifest_name] = {'bytes': manifest.stat().st_size,
                               'sha256': hashlib.sha256(manifest.read_bytes()).hexdigest()}
    seen = set()
    with tarfile.open(output, 'r|gz') as archive:
        for member in archive:
            if member.name not in expected or member.name in seen or not member.isfile():
                raise ValueError(f'Unexpected archive member: {member.name}')
            digest = hashlib.sha256()
            with archive.extractfile(member) as stream:
                for block in iter(lambda: stream.read(2**20), b''):
                    digest.update(block)
            record = expected[member.name]
            if member.size != record['bytes'] or digest.hexdigest() != record['sha256']:
                raise ValueError(f'Archive content mismatch: {member.name}')
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError('Archive member count mismatch')
    print(json.dumps({'archive': str(output), 'bytes': output.stat().st_size,
                      'files': len(files) + 1, 'all_member_hashes_verified': True}, ensure_ascii=False))


if __name__ == '__main__':
    main()
