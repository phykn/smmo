# SMMO

**SMMO** (Scattering-Matrix Method for Multilayer Optics) is a Python library for calculating exact optical properties (Transmission & Reflection) of multilayer structures, supporting both **coherent** interference and **incoherent** bulk media.

## Installation

```bash
pip install .
```

## Quick Start

```python
import numpy as np
from smmo import SMMO, make_config, make_layer

# 1. Light Configuration (Wavenumber cm⁻¹, Angle 0-90, Pol 's'/'p')
cfg = make_config(np.linspace(500, 2000, 100), 0, "s")

# 2. Layer Stack (n, k, thickness cm, coherence bool)
layers = [
    make_layer(np.full(100, 1.0), np.zeros(100), 0, False),     # Air
    make_layer(np.full(100, 1.5), np.zeros(100), 0.01, True),  # Thin Film
    make_layer(np.full(100, 1.0), np.zeros(100), 0, False)     # Air
]

# 3. Compute
res = SMMO(layers, cfg)()
print(f"T: {res['T'].mean():.4f}, R: {res['R'].mean():.4f}")
```

## Features
- **Numerical Stability**: Uses S-matrix to avoid exponential growth issues in TMM.
- **Mixed Coherence**: Handles thin films and thick substrates in a single stack.
- **Broadband**: Efficiently processes spectral data using NumPy.

## Citation

```bibtex
@article{lee2022machine,
  title={Machine learning analysis of broadband optical reflectivity of semiconductor thin film},
  author={Lee, Byeoungju and others},
  journal={Journal of the Korean Physical Society},
  year={2022}
}
```

## References
1. [Ko & Inkson, Phys. Rev. B 38.14 (1988)](https://doi.org/10.1103/PhysRevB.38.9945)
2. [Ko & Sambles, JOSA A 5.11 (1988)](https://doi.org/10.1364/JOSAA.5.001863)
