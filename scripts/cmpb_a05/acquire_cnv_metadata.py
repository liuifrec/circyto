"""Bounded public CNV table/metadata acquisition only; never raw sequence data."""
import datetime as dt
import json
from pathlib import Path
import sys
import urllib.request

from inspect_sources import digest


def main():
    root = Path(sys.argv[1])
    items = {}
    stem = "GSE278958_06_scRR-seq-DNA_IMR90_G1_summary_CNV_"
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278958/suppl/"
    for kind in ["states", "mappabilitynorm"]:
        name = stem + kind + "_bin50kb.txt.gz"
        items[name] = base + name
    items.update({
        "geo_series.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278958&targ=self&form=text&view=full",
        "geo_dna_sample.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8558875&targ=self&form=text&view=full",
        "Codes_Fig6_part1.txt": "https://raw.githubusercontent.com/mcbmieu/scRR-seq/0910c9a30ce1b7e8a56cb97ac07dcb63a585c626/MainFigures/Fig6/scripts/Codes_Fig6_part1.txt",
        "Codes_Fig6_part2.R": "https://raw.githubusercontent.com/mcbmieu/scRR-seq/0910c9a30ce1b7e8a56cb97ac07dcb63a585c626/MainFigures/Fig6/scripts/Codes_Fig6_part2.R",
    })
    records = []
    for name, url in items.items():
        path = root / "inputs" / name
        rec = {"name": name, "url": url, "accessed_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                data = response.read(3 * 1024**2 + 1)
                assert len(data) <= 3 * 1024**2, "public-file size cap exceeded"
                rec["response_url"] = response.url
            with path.open("xb") as f:
                f.write(data)
            rec.update(status="retrieved", bytes=len(data), sha256=digest(path))
        except Exception as exc:
            rec.update(status="unavailable", error=repr(exc))
        records.append(rec)
        print(json.dumps(rec), flush=True)
        (root / "records/public_acquisition.json").write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
