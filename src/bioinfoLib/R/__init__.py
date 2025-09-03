import os
import subprocess

from dotenv import load_dotenv

load_dotenv()
r_ld_path = (
    subprocess.check_output(["python", "-m", "rpy2.situation", "LD_LIBRARY_PATH"])
    .decode("utf-8")
    .strip()
)
# sometimes there are warning messages, so r_ld_path need to be cleaned
lines = r_ld_path.split("\n")
valid_paths = []
for line in lines:
    paths = line.split(":")
    for path in paths:
        if os.path.exists(path):
            valid_paths.append(path)
r_ld_path = ":".join(valid_paths)
existing = os.environ.get("LD_LIBRARY_PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{r_ld_path}:{existing}"

from . import utils
