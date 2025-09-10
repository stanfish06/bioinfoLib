import io
import re
from contextlib import redirect_stdout

import juliapkg

packages = {
    "Ripserer": ("aa79e827-bd0b-42a8-9f10-2b302677a641", None),
    "StatsBase": ("2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91", None),
    "Combinatorics": ("861a8166-3701-5b0c-9a16-15d98fcdc6aa", None),
    "Graphs": ("86223c79-3864-5bf0-83f7-82e725a168b6", None),
    "SimpleWeightedGraphs": ("47aef6b3-ad0c-573a-a1e2-d07658019622", None),
}


def get_installed_packages():
    f = io.StringIO()
    with redirect_stdout(f):
        juliapkg.status()
    output = f.getvalue()
    pkg_list = []
    uuid_list = []
    lines = output.strip().split("\n")
    pkg_section_start = False
    for ln in lines:
        ln = ln.strip()
        if ln:
            if pkg_section_start:
                if ":" in ln:
                    pkg_name, pkg_info = ln.split(":", 1)
                    pkg_name = pkg_name.strip()
                    pkg_list.append(pkg_name)
                    uuid = re.search(r"'uuid':\s*'([^']+)'", pkg_info)
                    uuid = uuid.group(1) if uuid else None
                    uuid_list.append(uuid)
            if ln == "Packages:":
                pkg_section_start = True
    return (pkg_list, uuid_list)


current_packages = get_installed_packages()
print("Julia package information:")
for package, (uuid, version) in packages.items():
    add_pkg = False
    if package not in current_packages[0]:
        add_pkg = True
    else:
        idx = current_packages[0].index(package)
        installed_uuit = current_packages[1][idx]
        if installed_uuit != uuid:
            print(
                f"{pkg_name} installed but from a different uuid. Will update its uuid"
            )
            add_pkg = True
    if add_pkg:
        try:
            if version:
                juliapkg.add(package, uuid=uuid, version=version)
            else:
                juliapkg.add(package, uuid=uuid)
            print(f"{package} was just added")
        except Exception as e:
            print(f"Error installing {package}: {e}")
            raise e
    else:
        print(f"{package}")

juliapkg.resolve()
from . import utils
