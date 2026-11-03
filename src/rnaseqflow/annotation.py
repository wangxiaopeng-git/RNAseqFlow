"""Extract transcript-to-gene pairs from GTF or GFF3 annotations.

No silent gene/transcript ID canonicalization; exact ID matching is intentional.
"""
from __future__ import annotations
import csv
import gzip
import re
from pathlib import Path
from urllib.parse import unquote

from .configuration import ConfigError

GTF_ATTR = re.compile(r'([A-Za-z_][\w.-]*)\s+"([^"]+)"')
TRANSCRIPT_TYPES = {"mRNA", "transcript", "primary_transcript", "lnc_RNA", "tRNA", "rRNA", "ncRNA"}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_gtf(attrs: str) -> dict[str, str]:
    return dict(GTF_ATTR.findall(attrs))


def parse_gff3(attrs: str) -> dict[str, str]:
    return dict((unquote(k.strip()), unquote(v.strip())) for x in attrs.split(";")
                if "=" in x for k, v in [x.split("=", 1)])


def extract_tx2gene(annotation: Path, output: Path, feature_format: str = "auto") -> int:
    annotation = Path(annotation)
    if not annotation.is_file():
        raise ConfigError(f"Annotation does not exist: {annotation}")
    if feature_format == "auto":
        stem = annotation.name.lower().removesuffix(".gz")
        if stem.endswith((".gtf", ".gff", ".gff3")):
            feature_format = "gtf" if stem.endswith(".gtf") else "gff3"
        else:
            raise ConfigError("Cannot auto-detect annotation format; use --format gtf or gff3")
    if feature_format not in ("gtf", "gff3"):
        raise ConfigError("Annotation format must be gtf or gff3")
    pairs = {}
    with open_text(annotation) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ConfigError(f"Annotation row {line_number} must have exactly nine GTF/GFF3 fields")
            kind, text = parts[2], parts[8]
            if feature_format == "gtf":
                attributes = parse_gtf(text)
                transcript = attributes.get("transcript_id")
                gene = attributes.get("gene_id")
            else:
                if kind not in TRANSCRIPT_TYPES:
                    continue
                attributes = parse_gff3(text)
                transcript = attributes.get("ID")
                gene = attributes.get("Parent")
                if gene and "," in gene:
                    raise ConfigError(f"GFF3 transcript {transcript} has multiple gene parents; resolve manually")
            if transcript and gene:
                if transcript in pairs and pairs[transcript] != gene:
                    raise ConfigError(f"Transcript {transcript} maps to conflicting genes")
                pairs[transcript] = gene
    if not pairs:
        raise ConfigError("No transcript-to-gene relationships found; check format/attributes")
    output = Path(output)
    if output.exists():
        raise ConfigError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id"])
        writer.writerows(sorted(pairs.items()))
    return len(pairs)
