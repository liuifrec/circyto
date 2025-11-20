import requests
import pandas as pd
from pathlib import Path

MANIFEST = Path("manifest.tsv")
OUT_TXT = Path("fastq_urls.txt")

# 1) read subset runs from manifest
mf = pd.read_csv(MANIFEST, sep="\t")
subset_runs = set(mf["cell_id"].astype(str))

print(f"Runs from manifest: {len(subset_runs)}")

# 2) download ENA filereport for the whole project (ERP104512)
ena_url = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession=ERP104512"
    "&result=read_run"
    "&fields=run_accession,fastq_ftp"
    "&format=tsv&download=true"
)

print("Fetching ENA filereport…")
resp = requests.get(ena_url)
resp.raise_for_status()
lines = resp.text.strip().splitlines()
header = lines[0].split("\t")
run_idx = header.index("run_accession")
ftp_idx = header.index("fastq_ftp")

urls = []
for line in lines[1:]:
    cols = line.split("\t")
    run = cols[run_idx]
    if run in subset_runs:
        fastq_ftp = cols[ftp_idx]
        if not fastq_ftp:
            continue
        # fastq_ftp can be "path1;path2"
        for rel in fastq_ftp.split(";"):
            rel = rel.strip()
            if not rel:
                continue
            urls.append("ftp://" + rel)

print(f"Collected {len(urls)} FASTQ URLs")
OUT_TXT.write_text("\n".join(urls) + "\n", encoding="utf-8")
print(f"Wrote {OUT_TXT.resolve()}")
