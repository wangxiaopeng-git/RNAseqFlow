from pathlib import Path


def test_installed_assets_match_canonical_sources():
    root=Path(__file__).resolve().parents[1]
    for path in list((root/'scripts').glob('*.R')) + list((root/'scripts').glob('*.py')) + [root/'workflow/Snakefile',root/'config/config.yaml']:
        assert path.read_bytes() == (root/'src/rnaseqflow/assets'/path.relative_to(root)).read_bytes(), path
