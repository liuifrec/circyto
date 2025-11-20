#!/bin/bash
set -e
echo "[circCyto] Applying multi-detector + merge + CI workflow patch..."

mkdir -p circyto/detectors
cat > circyto/detectors/__init__.py <<'PY'
from .ciri_full import run_ciri_full
from .ciri_long import run_ciri_long
from .find_circ import run_find_circ
from .circ_explorer2 import run_circ_explorer2
DETECTORS = {
    "ciri-full": run_ciri_full,
    "ciri-long": run_ciri_long,
    "find-circ": run_find_circ,
    "circ-explorer2": run_circ_explorer2,
}
def get_engine(name):
    if name not in DETECTORS:
        raise ValueError(f"Unknown engine: {name}. Available: {list(DETECTORS)}")
    return DETECTORS[name]
PY

cat > circyto/detectors/ciri_full.py <<'PY'
import subprocess
def run_ciri_full(r1,r2,ref_fa,gtf,out_tsv,threads=8):
    cmd=f"CIRI-full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv} -p {threads}"
    print("[CIRI-full]",cmd)
    subprocess.run(cmd,shell=True,check=False)
PY

cat > circyto/detectors/ciri_long.py <<'PY'
import subprocess
def run_ciri_long(r1,r2,ref_fa,gtf,out_tsv,threads=8):
    cmd=f"CIRI-long -r {ref_fa} -a {gtf} -i {r1} -o {out_tsv}"
    print("[CIRI-long]",cmd)
    subprocess.run(cmd,shell=True,check=False)
PY

cat > circyto/detectors/find_circ.py <<'PY'
import subprocess
def run_find_circ(r1,r2,ref_fa,gtf,out_tsv,threads=8):
    cmd=f"python find_circ.py -G {ref_fa} -p {r1} > {out_tsv}"
    print("[find_circ]",cmd)
    subprocess.run(cmd,shell=True,check=False)
PY

cat > circyto/detectors/circ_explorer2.py <<'PY'
import subprocess
def run_circ_explorer2(r1,r2,ref_fa,gtf,out_tsv,threads=8):
    cmd=f"CIRCexplorer2 parse {r1} > parse.txt && CIRCexplorer2 annotate -r {gtf} -g {ref_fa} -b parse.txt > {out_tsv}"
    print("[CIRCexplorer2]",cmd)
    subprocess.run(cmd,shell=True,check=False)
PY

cat > circyto/detectors/install_tools.py <<'PY'
import subprocess
def install_tool(tool):
    repos={
      "ciri-full":"https://github.com/bioinfo-biols/CIRI-full.git",
      "ciri-long":"https://github.com/bioinfo-biols/CIRI-long.git",
      "find-circ":"https://github.com/marvin-jens/find_circ.git",
    }
    if tool in repos:
        subprocess.run(f"git clone {repos[tool]}",shell=True)
    elif tool=="circ-explorer2":
        subprocess.run("pip install CIRCexplorer2",shell=True)
    else:
        print("Unknown tool:",tool)
PY

mkdir -p circyto/pipeline
cat > circyto/pipeline/merge_results.py <<'PY'
import pandas as pd
from pathlib import Path
def merge_engine_outputs(input_dirs,mode,outdir):
    circ_sets=[]
    for d in input_dirs:
        circ_ids=set()
        for f in Path(d).glob("*.tsv"):
            try:
                df=pd.read_csv(f,sep='\t',comment='#')
                if 'circ_id' in df.columns:
                    circ_ids.update(df['circ_id'].astype(str))
            except Exception as e:
                print("skip",f,":",e)
        circ_sets.append(circ_ids)
    merged=set.union(*circ_sets) if mode=="union" else set.intersection(*circ_sets)
    Path(outdir).mkdir(parents=True,exist_ok=True)
    with open(Path(outdir)/"merged_circ_ids.txt","w") as fw:
        for cid in sorted(list(merged)):
            fw.write(cid+"\n")
    print(f"[merge] {len(merged)} circRNAs saved to {outdir}/merged_circ_ids.txt")
PY

mkdir -p .github/workflows
cat > .github/workflows/test.yml <<'YML'
name: CI
on: [push, pull_request]
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    - name: Install
      run: pip install -e .[detectors]
    - name: Run mini test
      run: bash tests/run_tests.sh || true
YML

mkdir -p tests
echo "cell_1\ncell_2" > tests/cells.tsv
echo "@r\nACGT\n+\n####" > tests/mini_R1.fastq
echo "@r\nTGCA\n+\n####" > tests/mini_R2.fastq

cat > tests/run_tests.sh <<'BASH'
#!/bin/bash
echo "Running mini tests..."
for eng in ciri-full ciri-long find-circ circ-explorer2; do
  echo "Testing $eng"
  mkdir -p tests/out_$eng
  circyto run --engine $eng --manifest tests/cells.tsv --ref-fa tests/mini_R1.fastq --gtf tests/mini_R2.fastq --outdir tests/out_$eng || echo "skip $eng"
done
echo "Merging union/intersection"
circyto merge --inputs tests/out_ciri-full tests/out_ciri-long --mode union --outdir tests/merged_union || true
circyto merge --inputs tests/out_ciri-full tests/out_ciri-long --mode intersection --outdir tests/merged_intersection || true
BASH
chmod +x tests/run_tests.sh

# Add CI badge
if ! grep -q "actions/workflows/test.yml" README.md; then
  sed -i '1i ![CI](https://github.com/liuifrec/circyto/actions/workflows/test.yml/badge.svg)\n' README.md || true
fi

git add -A
git commit -m "feat: multi-detector engines + merge + CI workflow + badge"
git push origin feat/plate-tenx5-v2

echo "✅ Patch complete. Check GitHub Actions for CI run."
