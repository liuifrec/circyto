import subprocess


def run_circ_explorer2(r1, r2, ref_fa, gtf, out_tsv, threads=8):
    cmd = f"CIRCexplorer2 parse {r1} > parse.txt && CIRCexplorer2 annotate -r {gtf} -g {ref_fa} -b parse.txt > {out_tsv}"
    print("[CIRCexplorer2]", cmd)
    subprocess.run(cmd, shell=True, check=False)
