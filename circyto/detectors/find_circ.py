import subprocess
def run_find_circ(r1,r2,ref_fa,gtf,out_tsv,threads=8):
    cmd=f"python find_circ.py -G {ref_fa} -p {r1} > {out_tsv}"
    print("[find_circ]",cmd)
    subprocess.run(cmd,shell=True,check=False)
