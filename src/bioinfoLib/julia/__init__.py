import juliapkg

packages = {
    "Ripserer": ("aa79e827-bd0b-42a8-9f10-2b302677a641", None),
    "StatsBase": ("2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91", None),
    "Combinatorics": ("861a8166-3701-5b0c-9a16-15d98fcdc6aa", None),
    "Graphs": ("86223c79-3864-5bf0-83f7-82e725a168b6", None),
    "SimpleWeightedGraphs": ("47aef6b3-ad0c-573a-a1e2-d07658019622", None)
}

current_packages = juliapkg.status()
for package, (uuid, version) in packages.items():
    if package not in current_packages:
        try:
            if version:
                juliapkg.add(package, uuid=uuid, version=version)
            else:
                juliapkg.add(package, uuid=uuid)
        except Exception as e:
            print(f"Error installing {package}: {e}")
            raise e
juliapkg.resolve()
from . import utils
