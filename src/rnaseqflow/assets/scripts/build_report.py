"""Generate a static report index; R produces the statistics and figures."""
import argparse
import html
import json
from pathlib import Path
from urllib.parse import quote
from rnaseqflow.provenance import artifact_path


def render(summary: dict) -> str:
    title = html.escape(str(summary.get("project", "RNAseqFlow")))
    rows = []
    for sec in summary.get("sections", []):
        heading = html.escape(str(sec.get("title", "Section")))
        detail = html.escape(str(sec.get("detail", "")))
        for filename in sec.get("files", []):
            artifact_path(Path.cwd(), str(filename))
        images = ''.join(
            f'<figure><img loading="lazy" src="../{html.escape(quote(str(p), safe="/"), quote=True)}" alt="{html.escape(str(p), quote=True)}"><figcaption>{html.escape(str(p))}</figcaption></figure>'
            for p in sec.get('files', []) if str(p).endswith(('.svg', '.png')))
        files = "".join(
            f'<li><a href="../{html.escape(str(p), quote=True)}">{html.escape(str(p))}</a></li>'
            for p in sec.get("files", [])
        )
        rows.append(f"<section><h2>{heading}</h2><p>{detail}</p>{images}<ul>{files}</ul></section>")
    metrics = html.escape(f"Samples: {summary.get('samples', '?')} · Genes retained: {summary.get('genes_after_filter', '?')}")
    contrasts = ''.join(f"<li>{html.escape(str(c.get('name')))}: UP {int(c.get('up', 0))}, DOWN {int(c.get('down', 0))}, FDR significant {int(c.get('fdr_significant', 0))}</li>" for c in summary.get('contrasts', []))
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title} | RNAseqFlow</title>
<style>body{{font:16px/1.6 system-ui,sans-serif;max-width:1040px;margin:2rem auto;padding:0 1rem;color:#233}}
header{{border-bottom:3px solid #14816b;padding-bottom:1rem}}section{{padding:1rem;margin:1rem 0;background:#f6faf9;border-radius:12px}}
a{{color:#086b7e}}ul{{columns:2}}li{{break-inside:avoid;margin:.4rem 0}}
img{{max-width:100%;height:auto}}figure{{margin:1rem 0}}small{{color:#667}}@media(max-width:600px){{ul{{columns:1}}}}</style></head>
<body><header><h1>{title}</h1><p>RNAseqFlow v0.2.0 · level: {html.escape(str(summary.get('level','')))}</p>
<p>{metrics}</p><ul>{contrasts}</ul><small>Generated from analysis outputs. Statistical significance and biological validity require review.</small></header>
{''.join(rows)}<footer><small>Inputs, versions, and R session metadata are recorded alongside the results.</small></footer></body></html>'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    with open(args.summary, encoding="utf-8") as f:
        summary = json.load(f)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(summary), encoding="utf-8")


if __name__ == "__main__":
    main()
