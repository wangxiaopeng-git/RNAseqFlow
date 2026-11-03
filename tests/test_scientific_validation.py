
import pytest
import yaml

from rnaseqflow.cli import main
from rnaseqflow.configuration import ConfigError, _check_tsv_map, validate_config


@pytest.mark.parametrize('body', ['T1\tG1\nT1\tG2\n', 'T1\t\n', '\tG1\n', 'T1\tG1\textra\n'])
def test_ambiguous_or_empty_map_rejected(tmp_path, body):
    path = tmp_path / 'map.tsv'
    path.write_text('transcript_id\tgene_id\n' + body)
    with pytest.raises(ConfigError):
        _check_tsv_map(path)


def test_repeated_identical_map_is_allowed(tmp_path):
    path = tmp_path / 'map.tsv'
    path.write_text('transcript_id\tgene_id\nT1\tG1\nT1\tG1\n')
    _check_tsv_map(path)


@pytest.fixture
def project(tmp_path):
    dest = tmp_path / 'project'
    assert main(['init', '--dest', str(dest)]) == 0
    return dest / 'config/config.yaml'


@pytest.mark.parametrize('key,value', [
    ('log2fc', -1), ('log2fc', float('nan')), ('min_count', 1.5),
    ('min_samples', 100), ('top_heatmap_genes', 0), ('lfc_shrinkage', 'unknown'),
    ('covariate_types', {'absent': 'numeric'}),
    ('covariate_types', {'condition': 'numeric'}),
    ('covariate_types', {'condition': 'bogus'}),
    ('covariate_types', {'batch': 'numeric'}),
    ('design', '~ batch * condition'),
])
def test_invalid_analysis_rejected(project, key, value):
    cfg = yaml.safe_load(project.read_text())
    cfg['analysis'][key] = value
    project.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ConfigError):
        validate_config(project)


def test_contrast_types_normalized(project):
    cfg = yaml.safe_load(project.read_text())
    cfg['analysis']['lfc_shrinkage'] = 'normal'
    project.write_text(yaml.safe_dump(cfg))
    validated, _, _ = validate_config(project)
    assert validated['analysis']['covariate_types']['condition'] == 'categorical'


def test_empty_duplicate_sample_columns_rejected(tmp_path):
    from rnaseqflow.configuration import read_samples
    sheet = tmp_path / 's.tsv'
    sheet.write_text('sample_id\tcondition\tcondition\nS1\tA\tB\n')
    with pytest.raises(ConfigError):
        read_samples(sheet)
