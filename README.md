# build and install wheel locally
- uv build
- pip install dist/name.whl (this could be done in any virtual env with pip)
# use environment locally after syncing environment
- uv run python -m ipykernel install --name='kernel name' --user
