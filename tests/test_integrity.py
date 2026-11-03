import json
from pathlib import Path

import pytest
import yaml

from rnaseqflow.provenance import write_manifest, verify_manifest


def make_run(tmp_path):
    cfg = tmp_path/'config.yaml'
    source=tmp_path/'counts.tsv'; source.write_text('gene_id\tS1\nG1\t1\n')
    cfg.write_text(yaml.safe_dump({'project':{'input_mode':'counts','counts_matrix':str(source)},'analysis':{}}))
    out=tmp_path/'results/analysis'; out.mkdir(parents=True)
    (out/'qc.svg').write_text('<svg/>')
    summary=out/'summary.json'
    summary.write_text(json.dumps({'project':'test','sections':[{'title':'QC','files':['qc.svg']}]}))
    return summary,cfg


def test_manifest_detects_missing_and_corrupt_optional_plot(tmp_path):
    summary,cfg=make_run(tmp_path)
    write_manifest(summary,cfg)
    manifest=summary.with_name('provenance.json')
    assert verify_manifest(manifest)==[]
    data=json.loads(manifest.read_text())
    assert data['inputs'][0]['sha256']
    summary.with_name('qc.svg').write_text('changed')
    assert any('qc.svg' in error for error in verify_manifest(manifest))
    summary.with_name('qc.svg').unlink()
    assert any('qc.svg' in error for error in verify_manifest(manifest))


@pytest.mark.parametrize('path',['../outside.svg','/tmp/outside.svg','https://example.com/x.svg'])
def test_manifest_rejects_escaping_paths(tmp_path,path):
    summary,cfg=make_run(tmp_path)
    summary.write_text(json.dumps({'sections':[{'files':[path]}]}))
    with pytest.raises(ValueError):
        write_manifest(summary,cfg)


def test_bad_manifest_is_not_complete(tmp_path):
    p=tmp_path/'provenance.json'
    for text in ['{}','null','{"outputs": []}','not json']:
        p.write_text(text)
        assert verify_manifest(p)


def test_provenance_resolves_snapshot_inputs_from_project_root(tmp_path):
    from rnaseqflow.provenance import input_paths
    config=tmp_path/'results/run_config.yaml'
    config.parent.mkdir()
    assert input_paths({'project':{'input_mode':'counts','counts_matrix':'data/counts.tsv'}}, config)[1] == tmp_path/'data/counts.tsv'


def test_only_active_route_inputs_are_recorded(tmp_path):
    from rnaseqflow.provenance import input_paths
    cfg={'project':{'input_mode':'counts','counts_matrix':'counts.tsv','tx2gene':'unused.tsv','transcript_fasta':'unused.fa'}}
    paths=input_paths(cfg,tmp_path/'config.yaml')
    assert paths == [tmp_path/'config.yaml',tmp_path/'counts.tsv']


def test_report_embeds_figures_and_escapes_metadata():
    from rnaseqflow.assets.scripts.build_report import render
    output=render({'project':'<unsafe>','samples':7,'genes_after_filter':8148,
                   'sections':[{'title':'QC','files':['qc/pca.svg']}],
                   'contrasts':[{'name':'A_vs_B','up':114,'down':110}]})
    assert '<img' in output and '../qc/pca.svg' in output
    assert '8148' in output and '114' in output
    assert '<unsafe>' not in output


def test_report_rejects_path_traversal():
    from rnaseqflow.assets.scripts.build_report import render
    with pytest.raises(ValueError):
        render({'sections':[{'files':['../../secret.svg']}]})
