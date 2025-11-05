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
