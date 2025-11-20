import pandas as pd
from pathlib import Path


def merge_engine_outputs(input_dirs, mode, outdir):
    circ_sets = []
    for d in input_dirs:
        circ_ids = set()
        for f in Path(d).glob("*.tsv"):
            try:
                df = pd.read_csv(f, sep="\t", comment="#")
                if "circ_id" in df.columns:
                    circ_ids.update(df["circ_id"].astype(str))
            except Exception as e:
                print("skip", f, ":", e)
        circ_sets.append(circ_ids)
    merged = set.union(*circ_sets) if mode == "union" else set.intersection(*circ_sets)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    with open(Path(outdir) / "merged_circ_ids.txt", "w") as fw:
        for cid in sorted(list(merged)):
            fw.write(cid + "\n")
    print(f"[merge] {len(merged)} circRNAs saved to {outdir}/merged_circ_ids.txt")
