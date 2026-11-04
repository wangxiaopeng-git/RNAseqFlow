"""Record run provenance and verify every advertised analysis artifact."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import yaml


def artifact_path(root: Path, name: str) -> Path:
    """Resolve only local relative paths, including symlink containment checks."""
    p = PurePosixPath(name)
    if not name or p.is_absolute() or '..' in p.parts or ':' in name or '\\' in name:
        raise ValueError(f'Unsafe artifact path: {name}')
    resolved = (root / name).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f'Artifact escapes analysis directory: {name}')
    return resolved


def file_record(path: Path, name: str | None = None) -> dict:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': name or str(path.resolve()), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def input_paths(config: dict, config_path: Path) -> list[Path]:
    from .configuration import project_root
    root = project_root(config_path)
    if config_path.name == 'run_config.yaml' and config_path.parent.name == 'results':
        root = config_path.resolve().parent.parent
    def resolve(value):
        p = Path(value).expanduser()
        return p.resolve() if p.is_absolute() else (root / p).resolve()
    project = config['project']
    paths = [config_path.resolve()]
    keys = {'counts': ('samples', 'counts_matrix'), 'salmon': ('samples', 'tx2gene'),
            'fastq': ('samples', 'tx2gene', 'transcript_fasta')}[project['input_mode']]
    for key in keys:
        if project.get(key):
            paths.append(resolve(project[key]))
    if config.get('analysis', {}).get('level') != 'basic' and config.get('analysis', {}).get('enrichment', {}).get('enabled') and project.get('gene_sets_gmt'):
        paths.append(resolve(project['gene_sets_gmt']))
    if project.get('samples') and project.get('input_mode') in ('fastq', 'salmon'):
        with resolve(project['samples']).open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if project['input_mode'] == 'salmon':
                    paths.append(resolve(row['quant_dir']) / 'quant.sf')
                else:
                    paths.extend(resolve(row[k]) for k in ('fastq_1', 'fastq_2') if row.get(k))
                    paths.append(root / 'results/quant' / row['sample_id'] / 'quant.sf')
    return list(dict.fromkeys(paths))


def write_manifest(summary_path: Path, config_path: Path) -> None:
    summary_path, config_path = Path(summary_path), Path(config_path)
    root = summary_path.parent
    summary = json.loads(summary_path.read_text())
    config = yaml.safe_load(config_path.read_text())
    sections = summary.setdefault('sections', [])
    if not any('provenance.json' in s.get('files', []) for s in sections):
        sections.append({'title': 'File integrity and provenance',
                         'detail': 'SHA-256 input/output checksums, effective configuration, software and source versions.',
                         'files': ['provenance.json']})
    names = sorted({'summary.json'} | {f for s in sections for f in s.get('files', []) if f != 'provenance.json'})
    # Validate all paths and existence before updating the summary.
    for name in names:
        if not artifact_path(root, name).is_file():
            raise ValueError(f'Missing analysis artifact: {name}')
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
    versions = {'python': platform.python_version()}
    for package in ('rnaseqflow', 'PyYAML', 'snakemake'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = 'not installed'
    for tool in ('Rscript', 'fastp', 'fastqc', 'multiqc', 'salmon'):
        if shutil.which(tool):
            try:
                result = subprocess.run([tool, '--version'], capture_output=True, text=True, timeout=10)
                versions[tool] = (result.stdout + result.stderr).strip()[:1000]
            except (OSError, subprocess.TimeoutExpired):
                versions[tool] = 'version query unavailable'
    package_root = Path(__file__).parent
    code = [p for p in package_root.rglob('*') if p.is_file() and p.suffix in ('.R', '.py')]
    code += [package_root / 'assets/workflow/Snakefile']
    manifest = {'schema_version': 1, 'created_at': datetime.now(timezone.utc).isoformat(),
                'config': config, 'software': versions,
                'inputs': [file_record(p) for p in input_paths(config, config_path)],
                'code': [file_record(p) for p in sorted(set(code))],
                'outputs': [file_record(artifact_path(root, n), n) for n in names]}
    target = root / 'provenance.json'
    temporary = target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(target)


def verify_manifest(path: Path) -> list[str]:
    path = Path(path)
    try:
        manifest = json.loads(path.read_text())
        if manifest['schema_version'] != 1 or not manifest['outputs']:
            return ['Invalid or empty integrity manifest']
        errors = []
        names = set()
        for record in manifest['outputs']:
            name = record['path']
            names.add(name)
            target = artifact_path(path.parent, name)
            if not target.is_file() or file_record(target, name) != record:
                errors.append(f'Missing or changed artifact: {name}')
        summary = json.loads((path.parent / 'summary.json').read_text())
        required = {'summary.json'} | {f for s in summary['sections'] for f in s.get('files', []) if f != 'provenance.json'}
        if names != required:
            errors.append('Manifest does not cover all advertised outputs')
        return errors
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return ['Missing or invalid integrity manifest']


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    write_manifest(args.summary, args.config)


if __name__ == '__main__':
    main()
