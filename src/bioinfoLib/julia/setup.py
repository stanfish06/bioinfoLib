import juliapkg


def install_packages():
    packages = {
        "Ripserer": "aa79e827-bd0b-42a8-9f10-2b302677a641",
        "StatsBase": "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91",
        "Plots": "91a5bcdd-55d7-5caf-9e0b-520d859cae80",
    }
    juliapkg.require_julia("1.11")
    for package, uuid in packages.items():
        try:
            juliapkg.add(package, uuid=uuid)
        except Exception as e:
            print(f"Error installing {package}: {e}")
            raise e
    juliapkg.resolve()
