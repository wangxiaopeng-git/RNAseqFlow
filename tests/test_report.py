import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from build_report import render


def test_report_escapes_markup():
    output = render({'project':'<script>alert(1)</script>', 'level':'basic',
                     'sections':[{'title':'test', 'files':['a<b.tsv']} ]})
    assert '<script>' not in output
    assert '&lt;script&gt;' in output
    assert 'a&lt;b.tsv' in output
