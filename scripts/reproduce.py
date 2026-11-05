"""Rebuild the SL experiment in a separate directory, without using fitted models."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            value.update(block)
    return value.hexdigest()


def fixed_inputs(allow_modified=False):
    manifest = json.loads((ROOT / 'data/reference/input_manifest.json').read_text())
    records = []
    changed = []
    for expected in manifest['files']:
        path = ROOT / expected['path']
        if not path.resolve().is_relative_to(ROOT) or not path.is_file():
            raise FileNotFoundError(f'Missing fixed input: {expected["path"]}. Run git lfs pull after cloning.')
        actual = digest(path)
        if actual != expected['sha256']:
            changed.append(expected['path'])
        records.append({'path': expected['path'], 'sha256': actual, 'bytes': path.stat().st_size})
    if changed and not allow_modified:
        raise ValueError('Fixed inputs differ from the reference: ' + ', '.join(changed)
                         + '. Restore the reference inputs, or explicitly use --allow-modified-inputs for a new experiment.')
    return records, changed


def source_files():
    paths = []
    for folder in ['metagene', 'scripts']:
        paths.extend(p for p in (ROOT / folder).iterdir() if p.is_file() and p.suffix in {'.py', '.R'})
    paths.extend((ROOT / 'configs').glob('*.json'))
    return sorted(paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, help='New work directory; reference results are kept intact.')
    parser.add_argument('--devices', nargs='+', help='For example cuda:0 cuda:1, cuda:0, or cpu; automatic when omitted.')
    parser.add_argument('--config', type=Path, help='Experiment JSON to use instead of the reference configuration.')
    parser.add_argument('--check-only', action='store_true', help='Verify the nine fixed inputs without running models.')
    parser.add_argument('--prepare-only', action='store_true', help='Create a clean work directory without training.')
    parser.add_argument('--resume', action='store_true', help='Resume this work directory with identical code, inputs and configuration.')
    parser.add_argument('--allow-modified-inputs', action='store_true', help='Record changed input hashes for a new experiment.')
    args = parser.parse_args()
    inputs, changed = fixed_inputs(args.allow_modified_inputs)
    if args.check_only:
        print(json.dumps({'checked_inputs': len(inputs), 'bytes': sum(x['bytes'] for x in inputs),
                          'modified_inputs': changed, 'pretrained_weights_required': False}, ensure_ascii=False))
        return
    config = (args.config or ROOT / 'configs/experiment_sl.json').resolve()
    json.loads(config.read_text())
    sources = source_files()
    fingerprint = {'inputs': inputs, 'source': [{'path': p.relative_to(ROOT).as_posix(), 'sha256': digest(p)}
                                              for p in sources], 'experiment_config_sha256': digest(config)}
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    work = (args.output_dir or ROOT / 'reproductions' / stamp).resolve()
    if work == ROOT or ROOT.is_relative_to(work):
        raise ValueError('The reproduction directory must not replace the source repository or its parent.')
    launch = work / 'reproduction_launch.json'
    if args.resume:
        if not launch.is_file() or json.loads(launch.read_text())['fingerprint'] != fingerprint:
            raise ValueError('Resume requires the same source, fixed inputs and configuration. Use a new --output-dir.')
    else:
        if work.exists() and any(work.iterdir()):
            raise ValueError('Output directory is not empty. Use a new directory or --resume.')
        work.mkdir(parents=True, exist_ok=True)
        paths = set(sources)
        paths.update(ROOT / row['path'] for row in inputs)
        paths.add(ROOT / 'data/reference/input_manifest.json')
        paths.update(p for p in (ROOT / 'docs').iterdir() if p.is_file())
        paths.update(p for p in (ROOT / 'data').rglob('*.json') if p.stat().st_size < 300_000)
        paths.update(p for p in ROOT.glob('requirements*.txt'))
        paths.add(ROOT / 'README.md')
        for name in ['metabolic_adjudication.csv', 'chromatin_adjudication.csv']:
            path = ROOT / 'data/evidence' / name
            if path.exists():
                paths.add(path)
        for path in sorted(paths):
            target = work / path.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
        shutil.copy2(config, work / 'configs/experiment_sl.json')
        (work / 'results').mkdir(exist_ok=True)
        launch.write_text(json.dumps({'fingerprint': fingerprint, 'source_root': str(ROOT),
                                     'modified_inputs': changed, 'copied_pretrained_weights': False,
                                     'copied_training_runs': False}, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'work_directory': str(work), 'mode': 'prepared' if args.prepare_only else 'train_from_scratch',
                      'pretrained_weights_copied': False}, ensure_ascii=False), flush=True)
    if args.prepare_only:
        return
    command = [sys.executable, '-u', str(work / 'scripts/run_pipeline.py')]
    if args.devices:
        command += ['--devices', *args.devices]
    env = {**os.environ, 'MPLCONFIGDIR': str(work / '.matplotlib-cache')}
    with (work / 'reproduction.log').open('a') as log:
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(line, end='', flush=True)
        status = process.wait()
    if status:
        raise SystemExit(status)
    audit = json.loads((work / 'results/integrity_audit.json').read_text())
    if not audit['passed']:
        raise RuntimeError('Reproduction completed with failed integrity checks.')
    print(json.dumps({'completed': True, 'work_directory': str(work), 'audit_checks': audit['count'],
                      'report': str(work / 'results/MetaGene实验汇报.md')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
