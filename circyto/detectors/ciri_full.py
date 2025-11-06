import subprocess


def run_ciri_full(r1, r2, ref_fa, gtf, out_tsv, threads=8):
    cmd = f"CIRI-full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv} -p {threads}"
    print("[CIRI-full]", cmd)
    subprocess.run(cmd, shell=True, check=False)
