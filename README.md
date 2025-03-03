# build and install wheel locally
- uv build
- pip install dist/name.whl (this could be done in any virtual env with pip)
# use environment locally after syncing environment
- uv run python -m ipykernel install --name='kernel name' --user
# for linux
- install [cusparselt](https://developer.nvidia.com/cusparselt-downloads?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=24.04&target_type=deb_local)