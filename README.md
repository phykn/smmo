**SMMO** (**S**cattering-**M**atrix method for **M**ultilayer **O**ptics) is scattering matrix method [1-3] implementation for python.

### Installation

```python
pip install smmo
```

### Parameters

`config`

-   `wavenumber`: The wavenumber of refractive index (or absorption coefficient) of layers. The unit is inverse centimeter (cm<sup>-1</sup>).
-   `incidence`: Angle of incidence ($0\leqq\theta<90$). The unit is arc degree.
-   `polarization`: Polarization of the incidence light. `s` is s-polarization and `p` is p-polarization.

`layer`

-   `n`: Refractive index of the layer
-   `k`: Absorption coefficient of the layer
-   `thickness`: Thickness of the layer. The unit is centimeter (cm).
-   `coherence`: Coherence in the layer. `True` is coherence and `False` is incoherence.

### Example

```python
import numpy as np
from smmo import make_config, make_layer, SMMO

config = make_config(
    wavenumber=np.arange(0, 10000, step=1000),
    incidence=0,
    polarization="s"
)

layers = [
    make_layer(
        n=np.full(10, 1),
        k=np.full(10, 0),
        thickness=0,
        coherence=False
    ),
    make_layer(
        n=np.full(10, 1.5),
        k=np.full(10, 0),
        thickness=0.01,
        coherence=True
    ),
    make_layer(
        n=np.full(10, 2),
        k=np.full(10, 0),
        thickness=0.05,
        coherence=False
    ),
    make_layer(
        n=np.full(10, 1),
        k=np.full(10, 0),
        thickness=0,
        coherence=False
    )
]

output = SMMO(layers, config)()
```

```python
>>> print(output)
{'T': array([0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]),
 'R': array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2])}
```

### Citation

```python
@article{lee2022machine,
  title={Machine learning analysis of broadband optical reflectivity of semiconductor thin film},
  author={Lee, Byeoungju and Yu, Kwangnam and Jeon, Jiwon and Choi, EJ},
  journal={Journal of the Korean Physical Society},
  pages={1--5},
  year={2022},
  publisher={Springer}
}
```

### Reference

[[1](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.38.9945)] Ko, D. Yuk Kei, and J. C. Inkson., Physical Review B 38.14 9945 (1988)  
[[2](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-5-11-1863)] Ko, D. Yuk Kei, and J. R. Sambles., JOSA A 5.11 1863-1866 (1988)  
[[3](https://spie.org/Publications/Proceedings/Paper/10.1117/12.862566?SSO=1)] Dyakov, Sergey A., et al., International Conference on Micro-and Nano-Electronics 2009 (2010)
