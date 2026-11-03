from pathlib import Path
import shutil
from unittest.mock import patch

from rnaseqflow.cli import main


def test_job_limit_passed_to_real_snakemake_resources(tmp_path):
    project=tmp_path/'example'
    assert main(['init','--dest',str(project)])==0
    yaml_file=project/'config/config.yaml'
    with patch('rnaseqflow.cli.shutil.which', return_value='/usr/bin/snakemake'), \
         patch('rnaseqflow.cli.subprocess.call',return_value=0) as call:
        assert main(['run','-c',str(yaml_file),'--cores','18','--jobs','5',
                     '--mem-mb','50000','--dry-run'])==0
        argv=call.call_args.args[0]
        assert argv[argv.index('--cores')+1]=='18'
        assert 'concurrent_jobs=5' in argv
        assert 'mem_mb=50000' in argv
        assert '--dry-run' in argv
        assert call.call_args.kwargs['cwd']==project
    snapshot=project/'results/run_config.yaml'
    assert snapshot.is_file()
    before=snapshot.stat().st_mtime_ns
    with patch('rnaseqflow.cli.shutil.which', return_value='/usr/bin/snakemake'), \
         patch('rnaseqflow.cli.subprocess.call',return_value=0):
        assert main(['run','-c',str(yaml_file),'--cores','18','--jobs','5','--dry-run'])==0
    assert snapshot.stat().st_mtime_ns==before
