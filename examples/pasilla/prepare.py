"""Prepare a small real RNA-seq example, verifying fixed upstream SHA-256 hashes."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request

import yaml

COMMIT = '9f5623228d098720b28fe041604dc7a786e78f75'
BASE = f'https://raw.githubusercontent.com/bioconductor-source/pasilla/{COMMIT}/'
SOURCES = {
    'counts.tsv': ('inst/extdata/pasilla_gene_counts.tsv', 'ea0dafbfcc600559644cfe7dd5cc8de809d631eb64ba3089aaa25c2fa0954dad'),
    'sample_annotation.csv': ('inst/extdata/pasilla_sample_annotation.csv', '1844f871ca5d750f41685c73b87bd7dda291ec0d52304c1afb46a17ae88db745'),
    'DESCRIPTION': ('DESCRIPTION', '31ba017882bd1d244edd9dfd41915f79e88b42113fda112ad8ba6266b41b8b15'),
}


def prepare(dest: Path) -> None:
    if dest.exists():
        raise ValueError(f'Destination exists: {dest}; use a new directory')
    # Download and verify everything before creating the destination.
    data = {}
    for name, (relative, expected) in SOURCES.items():
        with urllib.request.urlopen(BASE + relative, timeout=60) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError(f'Upstream checksum mismatch: {name}')
        data[name] = content
    (dest/'data').mkdir(parents=True)
    (dest/'config').mkdir()
    for name, content in data.items():
        (dest/'data'/name).write_bytes(content)
    records=[{'file':name,'url':BASE+rel,'sha256':sha} for name,(rel,sha) in SOURCES.items()]
    (dest/'data/sources.json').write_text(json.dumps(records,indent=2)+'\n')
    with (dest/'data/samples.tsv').open('w', newline='') as stream:
        writer=csv.writer(stream,delimiter='\t')
        writer.writerow(['sample_id','condition','library_type'])
        for row in csv.DictReader(io.StringIO(data['sample_annotation.csv'].decode())):
            writer.writerow([row['file'].removesuffix('fb'),row['condition'],row['type'].replace('-','_')])
    cfg={'project':{'name':'Real pasilla RNAi','input_mode':'counts','samples':'data/samples.tsv','counts_matrix':'data/counts.tsv'},
         'analysis':{'level':'intermediate','design':'~ library_type + condition',
                     'covariate_types':{'library_type':'categorical','condition':'categorical'},
                     'contrasts':[{'name':'treated_vs_untreated','factor':'condition','numerator':'treated','denominator':'untreated'}],
                     'min_count':10,'min_samples':3,'fdr':.05,'log2fc':1.,'lfc_shrinkage':'none',
                     'enrichment':{'enabled':False},'wgcna':{'enabled':False}}}
    (dest/'config/config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    print(f'Prepared 7 real pasilla samples in {dest.resolve()}')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dest',type=Path,default=Path('pasilla_project'))
    prepare(parser.parse_args().dest)
