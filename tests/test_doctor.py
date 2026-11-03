from rnaseqflow.cli import main
from unittest.mock import patch


def test_doctor_reports_missing_tools(capsys):
    with patch('rnaseqflow.cli.shutil.which',return_value=None):
        assert main(['doctor']) == 1
    assert '"snakemake": null' in capsys.readouterr().out
