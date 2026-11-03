"""Opt-in actual R/Snakemake integration: RNASEQFLOW_INTEGRATION=1 pytest tests/integration."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

pytestmark = pytest.mark.skipif(os.environ.get('RNASEQFLOW_INTEGRATION') != '1', reason='requires full R environment')


def run(*args, cwd=None):
    result = subprocess.run([sys.executable, '-m', 'rnaseqflow.cli', *map(str,args)],
                            cwd=cwd, text=True,capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr


def project(path):
    run('init','--dest',path)
    cfg=path/'config/config.yaml'
    c=yaml.safe_load(cfg.read_text());c['analysis']['enrichment']['enabled']=False
    cfg.write_text(yaml.safe_dump(c))
    return cfg,c


def test_advanced_and_missing_optional_plot_regenerates(tmp_path):
    cfg,c=project(tmp_path/'advanced')
    c['analysis']['level']='advanced';cfg.write_text(yaml.safe_dump(c))
    run('run','-c',cfg,'--cores','2','--jobs','1','--mem-mb','24000')
    out=cfg.parent.parent/'results/analysis'
    assert (out/'network/status.txt').read_text().startswith('COMPLETED')
    run('verify','--results',out)
    (out/'qc/pca.svg').unlink()
    dry=run('run','-c',cfg,'--cores','2','--dry-run')
    assert 'regeneration' in dry and 'rule analyze:' in dry
    run('run','-c',cfg,'--cores','2','--jobs','1','--mem-mb','24000')
    assert (out/'qc/pca.svg').is_file()
    run('verify','--results',out)


def test_numeric_categories_and_shrinkage_preserve_inference(tmp_path):
    cfg,c=project(tmp_path/'shrink')
    sheet=cfg.parent.parent/'demo/samples.tsv'
    text=sheet.read_text().replace('control','0').replace('treated','1')
    sheet.write_text(text)
    ct=c['analysis']['contrasts'][0];ct['numerator']='1';ct['denominator']='0'
    cfg.write_text(yaml.safe_dump(c))
    run('run','-c',cfg,'--cores','2','--jobs','1','--mem-mb','24000')
    out=cfg.parent.parent/'results/analysis'
    table=out/'contrasts/treated_vs_control/deg_all.tsv'
    original={r['gene_id']:r for r in csv.DictReader(table.open(),delimiter='\t')}
    c['analysis']['lfc_shrinkage']='normal';cfg.write_text(yaml.safe_dump(c))
    run('run','-c',cfg,'--cores','2','--jobs','1','--mem-mb','24000')
    revised={r['gene_id']:r for r in csv.DictReader(table.open(),delimiter='\t')}
    for gene,row in original.items():
        for col in ['pvalue','padj','status','log2FoldChange']:
            assert revised[gene][col] == row[col]
    assert any(r['log2FoldChange_shrunk'] != r['log2FoldChange'] for r in revised.values())
    assert (out/'contrasts/treated_vs_control/ma_shrunk.svg').is_file()
    run('verify','--results',out)
