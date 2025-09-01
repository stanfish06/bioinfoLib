import os
import subprocess
from . import utils

from dotenv import load_dotenv

load_dotenv()
r_ld_path = (
    subprocess.check_output(["python", "-m", "rpy2.situation", "LD_LIBRARY_PATH"])
    .decode("utf-8")
    .strip()
)
existing = os.environ.get("LD_LIBRARY_PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{r_ld_path}:{existing}"
