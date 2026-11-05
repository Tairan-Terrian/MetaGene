"""Device selection and input fingerprints for reproducible SL runs."""
from pathlib import Path
import hashlib
import json
import os
import re
import time


def resolve_devices(devices=None):
    """Return validated execution devices; CLI values override the environment."""
    import torch

    if devices is None:
        devices = os.environ.get('METAGENE_DEVICES', '').strip()
    if isinstance(devices, str):
        devices = re.split(r'[\s,]+', devices.strip()) if devices.strip() else []
    else:
        devices = [part for item in (devices or [])
                   for part in re.split(r'[\s,]+', item.strip()) if part]
    if not devices:
        devices = ([f'cuda:{i}' for i in range(min(2, torch.cuda.device_count()))]
                   if torch.cuda.is_available() else ['cpu'])
    resolved = []
    for value in devices:
        device = torch.device(value)
        if device.type == 'cpu':
            value = 'cpu'
        elif device.type == 'cuda':
            index = 0 if device.index is None else device.index
            if not torch.cuda.is_available() or index >= torch.cuda.device_count():
                raise ValueError(f'CUDA device {index} is unavailable; select an available device or cpu.')
            value = f'cuda:{index}'
        else:
            raise ValueError(f'Unsupported device {value}; use cpu or cuda:<index>.')
        if value not in resolved:
            resolved.append(value)
    return resolved


def device_name(device):
    import torch

    selected = torch.device(device)
    return torch.cuda.get_device_name(selected) if selected.type == 'cuda' else 'CPU'


def candidate_annotations_path(root):
    """Prefer the published static input; accept the earlier server snapshot."""
    root = Path(root)
    current = root / 'data/reference/candidate_annotations.csv'
    legacy = root / 'results/sl/candidate_annotation_snapshot.csv'
    if current.is_file():
        return current
    if legacy.is_file():
        return legacy
    raise FileNotFoundError(
        'Missing data/reference/candidate_annotations.csv. Restore the fixed candidate '
        'input distributed with this experiment; prediction outputs are not candidate inputs.')


def run_fingerprint(root, config):
    """Freeze all current SL inputs and implementation files, excluding outputs."""
    root = Path(root)
    paths = [root / name for name in [
        'data/sl/ht29_labels.csv',
        'data/pairs_sources/HT29_gene_lfc.csv',
        'data/pairs_sources/HT29_annotations.csv',
        'data/annotations/symbol_mapping.csv',
        'data/processed/model_genes.csv',
        'data/evidence/integrated_features.npz',
        'data/evidence/graphs.npz',
        'data/evidence/external/external_pairs_model_domain.csv',
        'data/evidence/genes_annotation.csv',
        'configs/experiment_sl.json',
    ]]
    paths.append(candidate_annotations_path(root))
    paths.extend((root / 'metagene').glob('*.py'))
    paths.extend(root / 'scripts' / name for name in [
        'prepare_sl_labels.py', 'run_sl_experiments.py', 'run_pipeline.py',
        'deliver_sl.py', 'analyze_sl.py', 'audit_sl.py',
    ])
    files = []
    for path in sorted(set(paths)):
        if not path.is_file():
            raise FileNotFoundError(f'Missing reproduction input or implementation: {path}')
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(2**20), b''):
                digest.update(block)
        files.append({'path': path.relative_to(root).as_posix(), 'sha256': digest.hexdigest()})
    return {'schema_version': 2, 'files': files, 'config': config}


def _write_json_atomic(path, value):
    temporary = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    try:
        temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def ensure_run_manifest(root, config, out):
    """Reject stale run artifacts, while permitting a clone of reference summaries."""
    root, out = Path(root), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    expected = run_fingerprint(root, config)
    path = out / 'input_manifest.json'
    has_run_artifacts = (out / 'runs').is_dir() and any(
        item.is_file() for item in (out / 'runs').rglob('*'))
    previous = json.loads(path.read_text()) if path.is_file() else None
    same = previous is not None and all(previous.get(key) == expected[key]
                                        for key in ['schema_version', 'files', 'config'])
    if has_run_artifacts and not same:
        raise ValueError(
            'Existing SL run artifacts do not match the current inputs, code or configuration. '
            'Use a new reproduction work directory; existing results have not been changed.')
    if not same:
        manifest = {**expected, 'started_unix': time.time()}
        _write_json_atomic(path, manifest)
        _write_json_atomic(out / 'frozen_config.json', config)
        return manifest
    return previous
