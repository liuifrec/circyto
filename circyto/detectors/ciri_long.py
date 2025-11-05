import subprocess
def run_ciri_long(r1,r2,ref_fa,gtf,out_tsv,threads=8):
    cmd=f"CIRI-long -r {ref_fa} -a {gtf} -i {r1} -o {out_tsv}"
    print("[CIRI-long]",cmd)
    subprocess.run(cmd,shell=True,check=False)
