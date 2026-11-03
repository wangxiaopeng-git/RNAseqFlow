"""Configuration parsing and preflight validation.

All file paths are relative to the directory containing config/config.yaml's
parent directory (the project root). Absolute paths are also accepted.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import yaml

SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
LEVELS = ("basic", "intermediate", "advanced")
MODES = ("counts", "salmon", "fastq")


class ConfigError(ValueError):
    pass


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        obj = yaml.safe_load(fh)
    if not isinstance(obj, dict):
        raise ConfigError("Configuration must be a YAML mapping")
    return obj


def project_root(config_path: str | Path) -> Path:
    path = Path(config_path).resolve()
    if path.parent.name == "config":
        return path.parent.parent
    return path.parent


def resolve(root: Path, path: str | None) -> Path:
    if not path or not isinstance(path, str):
        raise ConfigError("Missing file path in configuration")
    p = Path(path).expanduser()
    return p.resolve() if p.is_absolute() else (root / p).resolve()


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise ConfigError("Sample sheet requires 'sample_id' header")
        if any(not x for x in reader.fieldnames):
            raise ConfigError("Sample sheet contains an empty column name")
        rows = list(reader)
    if not rows:
        raise ConfigError("Sample sheet has no samples")
    ids = [x.get("sample_id", "").strip() for x in rows]
    if len(set(ids)) != len(ids):
        raise ConfigError("Duplicate sample_id values")
    for ident in ids:
        if not SAFE_SAMPLE.fullmatch(ident) or ident in (".", ".."):
            raise ConfigError(f"Invalid sample_id {ident!r}; use letters, digits, . _ -")
    for row in rows:
        for key, val in row.items():
            if key is None:
                raise ConfigError("Extra fields in a sample sheet row")
            row[key] = (val or "").strip()
    return rows


def _check_counts(path: Path, sample_ids: list[str]) -> None:
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "gene_id":
            raise ConfigError("Count matrix must start with 'gene_id' column")
        if len(header) != len(set(header)):
            raise ConfigError("Count matrix has duplicate column names")
        missing = set(sample_ids) - set(header[1:])
        extra = set(header[1:]) - set(sample_ids)
        if missing or extra:
            raise ConfigError(f"Count sample columns differ from sheet: missing={sorted(missing)} extra={sorted(extra)}")
        genes = set()
        n = 0
        for lineno, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ConfigError(f"Count matrix row {lineno} has wrong column count")
            if not row[0] or row[0] in genes:
                raise ConfigError(f"Empty or duplicated gene_id at row {lineno}")
            genes.add(row[0])
            for value in row[1:]:
                try:
                    number = int(value)
                except (ValueError, TypeError):
                    raise ConfigError(f"Non-integer count at row {lineno}: {value!r}") from None
                if number < 0:
                    raise ConfigError(f"Negative count at row {lineno}")
            n += 1
        if n < 2:
            raise ConfigError("Count matrix needs at least two genes")


def _check_tsv_map(path: Path) -> None:
    with path.open(encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header or header[:2] != ["transcript_id", "gene_id"]:
            raise ConfigError("tx2gene must have transcript_id and gene_id as its first two columns")
        if not next(reader, None):
            raise ConfigError("tx2gene is empty")


def _design_columns(design: str) -> set[str]:
    if not isinstance(design, str) or not design.strip().startswith("~"):
        raise ConfigError("analysis.design must be an R formula starting with ~")
    # Standard formula terms and interactions only; intentionally no R expressions.
    body = design.strip()[1:].strip()
    if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*(?:\s*[+*: ]\s*[A-Za-z_][A-Za-z_0-9]*)*(?:\s*[+]\s*[01])?|1", body):
        raise ConfigError("Use a simple design formula, e.g. '~ batch + condition' or '~ sex * condition'")
    return set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", body))


def validate_config(config_path: str | Path, level_override: str | None = None) -> tuple[dict, Path, list[dict]]:
    config_path = Path(config_path).resolve()
    cfg = load_config(config_path)
    root = project_root(config_path)
    project = cfg.get("project", {})
    analysis = cfg.get("analysis", {})
    if not isinstance(project, dict) or not isinstance(analysis, dict):
        raise ConfigError("project and analysis must be YAML mappings")
    mode = project.get("input_mode")
    if mode not in MODES:
        raise ConfigError(f"input_mode must be one of {MODES}")
    level = level_override or analysis.get("level", "basic")
    if level not in LEVELS:
        raise ConfigError(f"level must be one of {LEVELS}")
    analysis["level"] = level
    cfg["analysis"] = analysis
    sheet_path = resolve(root, project.get("samples"))
    if not sheet_path.is_file():
        raise ConfigError(f"Sample sheet does not exist: {sheet_path}")
    rows = read_samples(sheet_path)
    if mode == "counts":
        counts = resolve(root, project.get("counts_matrix"))
        if not counts.is_file():
            raise ConfigError(f"Count matrix does not exist: {counts}")
        _check_counts(counts, [row["sample_id"] for row in rows])
    elif mode == "salmon":
        for row in rows:
            qp = row.get("quant_dir")
            if not qp or not (resolve(root, qp) / "quant.sf").is_file():
                raise ConfigError(f"Missing quant.sf for {row['sample_id']}: {qp}")
        map_path = resolve(root, project.get("tx2gene"))
        if not map_path.is_file():
            raise ConfigError(f"Missing tx2gene: {map_path}")
        _check_tsv_map(map_path)
    else:
        transcript = resolve(root, project.get("transcript_fasta"))
        txmap = resolve(root, project.get("tx2gene"))
        if not transcript.is_file() or not txmap.is_file():
            raise ConfigError("fastq mode requires existing transcript_fasta and tx2gene")
        _check_tsv_map(txmap)
        is_paired = None
        for row in rows:
            if not row.get("fastq_1"):
                raise ConfigError(f"Missing fastq_1 for {row['sample_id']}")
            for key in ("fastq_1", "fastq_2"):
                if row.get(key) and not resolve(root, row[key]).is_file():
                    raise ConfigError(f"Missing {key} for {row['sample_id']}: {row[key]}")
            if is_paired is None:
                is_paired = bool(row.get("fastq_2"))
            elif is_paired != bool(row.get("fastq_2")):
                raise ConfigError("Mixing single-end and paired-end FASTQs is not supported in one project")
    if level in ("intermediate", "advanced"):
        cols = set(rows[0])
        design_cols = _design_columns(analysis.get("design", ""))
        if not design_cols <= cols:
            raise ConfigError(f"Design variables missing from samples.tsv: {sorted(design_cols - cols)}")
        for var in design_cols:
            if any(not row.get(var) for row in rows):
                raise ConfigError(f"Missing or empty design covariate: {var}")
        contrasts = analysis.get("contrasts", [])
        if not isinstance(contrasts, list) or not contrasts:
            raise ConfigError("At least one contrast is required at intermediate/advanced level")
        seen = set()
        for c in contrasts:
            if not isinstance(c, dict) or not all(c.get(k) for k in ("name", "factor", "numerator", "denominator")):
                raise ConfigError("Each contrast needs name, factor, numerator, denominator")
            if not SAFE_SAMPLE.fullmatch(str(c["name"])) or c["name"] in seen:
                raise ConfigError("Contrast names must be unique and safe for filenames")
            seen.add(c["name"])
            if c["factor"] not in design_cols:
                raise ConfigError(f"Contrast factor {c['factor']} must appear in design")
            if c["numerator"] == c["denominator"]:
                raise ConfigError("Contrast numerator and denominator must differ")
            for value in (c["numerator"], c["denominator"]):
                replicates = sum(row[c["factor"]] == str(value) for row in rows)
                if replicates < 2:
                    raise ConfigError(f"Contrast {c['name']} has fewer than 2 samples in group {value}")
    for key in ("min_count", "min_samples"):
        if analysis.get(key, 1) < 1:
            raise ConfigError(f"{key} must be positive")
    fdr = float(analysis.get("fdr", 0.05))
    if not 0 < fdr < 1:
        raise ConfigError("fdr must be between 0 and 1")
    gmt = project.get("gene_sets_gmt")
    if level != "basic" and analysis.get("enrichment", {}).get("enabled", False):
        if not gmt or not resolve(root, gmt).is_file():
            raise ConfigError("Enrichment enabled but project.gene_sets_gmt is absent")
    if level == "advanced" and analysis.get("wgcna", {}).get("enabled", False):
        traits = analysis.get("wgcna", {}).get("traits", [])
        if not traits or not set(traits) <= set(rows[0]):
            raise ConfigError("wgcna.traits requires sample metadata columns")
    return cfg, root, rows
