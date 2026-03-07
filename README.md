# SMMO (Scattering-Matrix Method for Multilayer Optics)

**SMMO** is a Python implementation of the scattering-matrix method [1-3] designed to calculate the exact optical properties (Transmission and Reflection) of multilayer thin-film structures. 

### Background & Purpose
In optics, predicting how light propagates through a stack of different materials is crucial for designing anti-reflective coatings, optical filters, and analyzing semiconductor thin films. While traditional transfer-matrix methods are common, they often suffer from numerical instability at large thicknesses or high absorption conditions in the layers. 

SMMO employs the mathematically robust **Scattering-Matrix Method**, which resolves these mathematical instabilities. Furthermore, real-world optical measurements often involve a mix of thin microscopic layers (where interference and phase coherence matter) and thick macroscopic substrates (where incoherence dominates). SMMO explicitly models both **coherent** (phase-preserving) and **incoherent** layers seamlessly in a single calculation stack, making it a reliable and ideal toolkit for analyzing broadband optical reflectivity in experiments.

### Installation

To install SMMO from the local source directory, navigate to the project folder and run:

```bash
pip install .
```

### Usage

The `smmo` module operates on two main components to run its calculations: a `config` dictionary that dictates incident light properties, and a list of `layer` dictionaries that outline the physical structure.

#### Configuration (`make_config`)
Defines the parameters of the incident external light source.
- `wavenumber`: Mathematical inverse of the wavelength (cm<sup>-1</sup>). Provided as a numpy array.
- `incidence`: Angle of incidence $(0\leqq\theta<90)$ measured in arc degrees.
- `polarization`: The polarization state of the light. Use `"s"` for s-polarization and `"p"` for p-polarization.

#### Layers (`make_layer`)
Defines the physical and optical characteristics of each structural layer, ordered from top (incident side) to bottom.
- `n`: Real part of the refractive index. Must be a numpy array equal in length to the wavenumbers.
- `k`: Imaginary part of the refractive index (absorption coefficient). Must match the length of `n`.
- `thickness`: Physical thickness of the target layer in centimeters (cm).
- `coherence`: Boolean flag (`True` / `False`). Set to `True` if the layer is thin enough to trigger optical interference (coherent layer), and `False` for thick bulk media like substrates (incoherent layer).

### Example

The following script calculates the broadband optical reflection and transmission of a 4-layer structure involving a thin coherent film sandwiched between incoherent air and substrate layers.

```python
import numpy as np
from smmo import SMMO, make_config, make_layer

# 1. Set up the incidence condition
config = make_config(
    wavenumber=np.arange(0, 10000, step=1000),
    incidence=0,
    polarization="s"
)

# 2. Build the structural stack (from top light entry to bottom exit)
layers = [
    make_layer(n=np.full(10, 1.0), k=np.full(10, 0.0), thickness=0, coherence=False),      # Air
    make_layer(n=np.full(10, 1.5), k=np.full(10, 0.0), thickness=0.01, coherence=True),    # Coherent Film
    make_layer(n=np.full(10, 2.0), k=np.full(10, 0.0), thickness=0.05, coherence=False),   # Incoherent Substrate
    make_layer(n=np.full(10, 1.0), k=np.full(10, 0.0), thickness=0, coherence=False)       # Air
]

# 3. Compute Transmission (T) and Reflection (R)
output = SMMO(layers, config)()
print("Transmission:", output["T"])
print("Reflection:", output["R"])
```

### Citation

If you use this project in your research, please cite:
```bibtex
@article{lee2022machine,
  title={Machine learning analysis of broadband optical reflectivity of semiconductor thin film},
  author={Lee, Byeoungju and Yu, Kwangnam and Jeon, Jiwon and Choi, EJ},
  journal={Journal of the Korean Physical Society},
  pages={1--5},
  year={2022},
  publisher={Springer}
}
```

### References

1. [Ko, D. Yuk Kei, and J. C. Inkson., Physical Review B 38.14 9945 (1988)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.38.9945)
2. [Ko, D. Yuk Kei, and J. R. Sambles., JOSA A 5.11 1863-1866 (1988)](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-5-11-1863)
3. [Dyakov, Sergey A., et al., International Conference on Micro-and Nano-Electronics 2009 (2010)](https://spie.org/Publications/Proceedings/Paper/10.1117/12.862566?SSO=1)
