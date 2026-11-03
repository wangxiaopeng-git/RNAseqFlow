from pathlib import Path
import csv
import shutil

import pytest

from rnaseqflow.configuration import ConfigError, validate_config, _design_columns
from rnaseqflow.cli import main

ROOT = Path(__file__).resolve().parents[1]


def test_demo_valid():
    cfg, root, rows = validate_config(ROOT / "config/config.yaml")
    assert cfg["project"]["input_mode"] == "counts"
    assert len(rows) == 24
    assert root == ROOT


def test_design():
    assert _design_columns("~ batch + condition") == {"batch", "condition"}
    assert _design_columns("~ sex * condition") == {"sex", "condition"}
    with pytest.raises(ConfigError):
        _design_columns("~ system('rm -rf /')")


def test_init_and_validate(tmp_path):
    dest = tmp_path / "newproj"
    assert main(["init", "--dest", str(dest)]) == 0
    assert main(["validate", "-c", str(dest / "config/config.yaml")]) == 0
    assert main(["init", "--dest", str(dest)]) == 2


def test_invalid_counts(tmp_path):
    dest = tmp_path / "project"
    assert main(["init", "--dest", str(dest)]) == 0
    path = dest / "demo/gene_counts.tsv"
    text = path.read_text()
    assert "gene_id" in text
    path.write_text(text.replace("gene_id", "wrong_id", 1))
    with pytest.raises(ConfigError, match="gene_id"):
        validate_config(dest / "config/config.yaml")


def test_missing_replicates(tmp_path):
    dest = tmp_path / "project"
    assert main(["init", "--dest", str(dest)]) == 0
    import yaml
    c = dest / "config/config.yaml"
    cfg = yaml.safe_load(c.read_text())
    cfg["analysis"]["contrasts"][0]["numerator"] = "absent"
    c.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ConfigError, match="fewer than 2"):
        validate_config(c)


def test_salmon_missing_map(tmp_path):
    dest = tmp_path / "project"
    assert main(["init", "--dest", str(dest)]) == 0
    import yaml
    c = dest / "config/config.yaml"
    cfg = yaml.safe_load(c.read_text())
    cfg["project"]["input_mode"] = "salmon"
    c.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ConfigError, match="quant.sf"):
        validate_config(c)
