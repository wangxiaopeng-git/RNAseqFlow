import pytest
from rnaseqflow.annotation import extract_tx2gene
from rnaseqflow.configuration import ConfigError


def test_gtf(tmp_path):
    a = tmp_path / "test.gtf"
    a.write_text('chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                 'chr1\tsrc\texon\t11\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                 'chr1\tsrc\texon\t1\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n')
    out = tmp_path / 'pairs.tsv'
    assert extract_tx2gene(a,out) == 2
    assert out.read_text().splitlines() == ['transcript_id\tgene_id','T1\tG1','T2\tG1']
    with pytest.raises(ConfigError,match='already exists'):
        extract_tx2gene(a,out)


def test_gff3(tmp_path):
    a = tmp_path / "genes.gff3"
    a.write_text('##gff-version 3\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=G1\n'
                 'chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=T1;Parent=G1\n')
    out = tmp_path / 'pairs.tsv'
    assert extract_tx2gene(a,out) == 1
    assert 'T1\tG1' in out.read_text()


def test_conflicting_mapping(tmp_path):
    a = tmp_path / "wrong.gtf"
    a.write_text('chr1\tx\texon\t1\t5\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                 'chr1\tx\texon\t6\t10\t.\t+\t.\tgene_id "G2"; transcript_id "T1";\n')
    with pytest.raises(ConfigError,match="conflicting"):
        extract_tx2gene(a,tmp_path/'out.tsv')
