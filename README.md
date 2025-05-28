# bioinfoLib

A comprehensive Python library for bioinformatics analysis, focusing on single-cell RNA sequencing, image analysis, and Gaussian Process modeling.

## Features

- **Single-cell RNA Sequencing Analysis**
  - Comprehensive tools for scRNA-seq data processing and analysis
  - Integration with popular frameworks like Scanpy and CellRank

- **Image Analysis**
  - Advanced image processing and analysis capabilities
  - Integration with scikit-image and other imaging libraries

- **Gaussian Process Modeling**
  - Custom GP implementations for biological data modeling
  - Efficient Cython-based simulator

## Installation

### Prerequisites

- Python >= 3.12
- CUDA support for GPU acceleration (recommended)

### Local Installation

1. Clone the repository:
```bash
git clone https://github.com/stanfish06/bioinfoLib.git
cd bioinfoLib
```

2. Create and activate a virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. Install the package:
```bash
uv build
pip install dist/bioinfoLib-*.whl
```

### Jupyter Integration

To use bioinfoLib in Jupyter notebooks:

```bash
uv run python -m ipykernel install --name='bioinfoLib' --user
```

### Linux-specific Dependencies

On Linux systems, you may need to install additional CUDA dependencies:
```bash
# Install cusparselt (for Linux)
# Download from: https://developer.nvidia.com/cusparselt-downloads
# Install cuDNN
# Download from https://developer.nvidia.com/cudnn-downloads
# Install nccl
# Download from https://developer.nvidia.com/nccl/nccl-download

# on hpc, if you cannot install those manually, you can use conda to install them and set env variable to the conda lib folder
conda install nvidia::libcusparse conda-force::nccl nvidia::cuda-toolkit
# then
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

## Usage

```python
import bioinfoLib as binf

# Single-cell RNA sequencing analysis
from bioinfoLib import scRNAseq as scr

# Image analysis
from bioinfoLib import imageAnalysis as img

# Gaussian Process modeling
from bioinfoLib import GP

# Install julia dependencies.
bio.jl.setup.install_packages()
```

## Documentation

For detailed documentation and examples, please visit our [documentation page](https://github.com/stanfish06/bioinfoLib/wiki).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

[Add your license information here]

## Contact

[Add your contact information here]
