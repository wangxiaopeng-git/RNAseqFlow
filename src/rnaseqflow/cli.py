"""RNAseqFlow command-line interface."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path

import yaml

from .configuration import ConfigError, validate_config
from .annotation import extract_tx2gene
from .provenance import verify_manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="rnaseqflow", description="Modular bulk RNA-seq workflows")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("doctor", help="Check tools and R packages in the current environment")
    verify = sub.add_parser("verify", help="Verify hashes of all advertised analysis outputs")
    verify.add_argument("--results", type=Path, default=Path("results/analysis"))
    init = sub.add_parser("init", help="Create a runnable example project")
    init.add_argument("--dest", default="rnaseq_project")
    val = sub.add_parser("validate", help="Validate config and input data")
    val.add_argument("-c", "--config", default="config/config.yaml")
    run = sub.add_parser("run", help="Run or dry-run with Snakemake")
    run.add_argument("-c", "--config", default="config/config.yaml")
    run.add_argument("--level", choices=("basic", "intermediate", "advanced"))
    run.add_argument("--cores", type=int, default=8)
    run.add_argument("--jobs", type=int, default=4, help="Maximum simultaneous Snakemake jobs")
    run.add_argument("--mem-mb", type=int, default=32000, help="Global memory capacity in MB")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--keep-going", action="store_true")
    annotate = sub.add_parser("tx2gene", help="Build transcript-to-gene mapping from GTF/GFF3")
    annotate.add_argument("-a", "--annotation", required=True)
    annotate.add_argument("-o", "--output", required=True)
    annotate.add_argument("--format", choices=("auto","gtf","gff3"), default="auto")
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            errors = verify_manifest(args.results / 'provenance.json')
            print(json.dumps({'valid': not errors, 'errors': errors}, indent=2))
            return 1 if errors else 0
        if args.command == "doctor":
            binaries = ["snakemake", "Rscript", "fastp", "fastqc", "multiqc", "salmon"]
            present = {name: shutil.which(name) for name in binaries}
            if present["Rscript"]:
                rcode = 'p<-c("DESeq2","tximport","fgsea","WGCNA","yaml","ggplot2","pheatmap","jsonlite","svglite");cat(paste(p,vapply(p,requireNamespace,logical(1),quietly=TRUE),sep="=",collapse="\n"))'
                check = subprocess.run([present["Rscript"], "-e", rcode], capture_output=True, text=True)
                present["R_packages"] = check.stdout.strip() if check.returncode == 0 else f"ERROR: {check.stderr.strip()}"
            print(json.dumps(present, indent=2))
            r_ok = (present.get("R_packages") is not None and
                    all(f"{name}=TRUE" in present["R_packages"]
                        for name in ("DESeq2", "tximport", "fgsea", "WGCNA", "yaml", "ggplot2", "pheatmap", "jsonlite", "svglite")))
            return 0 if all(present[name] for name in binaries) and r_ok else 1
        if args.command == "tx2gene":
            count = extract_tx2gene(Path(args.annotation), Path(args.output), args.format)
            print(f"Wrote {count} transcript mappings to {args.output}")
            return 0
        if args.command == "init":
            dest = Path(args.dest).resolve()
            if dest.exists():
                raise ConfigError(f"Destination exists: {dest}. Choose a new directory.")
            assets = resources.files("rnaseqflow").joinpath("assets")
            dest.mkdir(parents=True)
            for folder in ("config", "demo"):
                source = assets.joinpath(folder)
                shutil.copytree(str(source), dest / folder)
            print(f"Created {dest}. Try: rnaseqflow validate -c {dest / 'config/config.yaml'}")
            return 0
        config = Path(args.config).resolve()
        cfg, root, samples = validate_config(config, getattr(args, "level", None))
        if args.command == "validate":
            print(json.dumps({"valid": True, "root": str(root), "mode": cfg["project"]["input_mode"],
                              "samples": len(samples), "level": cfg["analysis"]["level"]}, indent=2))
            return 0
        if min(args.cores, args.jobs, args.mem_mb) <= 0:
            raise ConfigError("--cores, --jobs, --mem-mb must be positive")
        exe = shutil.which("snakemake")
        if not exe:
            raise ConfigError("snakemake command not found. Install envs/python.yaml or pip install 'rnaseqflow[workflow]'")
        assets = resources.files("rnaseqflow").joinpath("assets")
        snakefile = Path(str(assets.joinpath("workflow/Snakefile")))
        # Snapshot config in user project for reproducibility; original is never edited.
        root.joinpath("results").mkdir(exist_ok=True)
        snapshot = root / "results" / "run_config.yaml"
        content = yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False)
        if not snapshot.is_file() or snapshot.read_text(encoding="utf-8") != content:
            snapshot.write_text(content, encoding="utf-8")
        cmd = [exe, "--snakefile", str(snakefile), "--configfile", str(snapshot),
               "--cores", str(args.cores), "--resources",
               f"concurrent_jobs={args.jobs}", f"mem_mb={args.mem_mb}",
               "--rerun-incomplete", "--printshellcmds"]
        analysis_dir = root / 'results/analysis'
        if (analysis_dir / 'summary.json').exists():
            errors = verify_manifest(analysis_dir / 'provenance.json')
            if errors:
                print('Analysis integrity check requires regeneration: ' + '; '.join(errors), flush=True)
                cmd.extend(['--forcerun', 'analyze'])
        if args.dry_run:
            cmd.append("--dry-run")
        if args.keep_going:
            cmd.append("--keep-going")
        print("Executing:", " ".join(cmd), flush=True)
        return subprocess.call(cmd, cwd=root)
    except (ConfigError, OSError, yaml.YAMLError) as exc:
        print(f"RNAseqFlow error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
