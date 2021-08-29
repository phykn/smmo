### Scattering-Matrix method for multilayer optics
Transfer-matrix method(TMM) is a simple, accurate and fast method used to understand the propagation of electromagnetic waves in a multilayer structure. However, TMM is numerically unstable when a thick layer in the structure has high absorption because exponential parts cause overflow. To overcome this problem, scattering-matrix method(SMM) was presented ([PRB-1988](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.38.9945), [JOSAA-1988](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-5-11-1863), [SPIE-2010](https://spie.org/Publications/Proceedings/Paper/10.1117/12.862566?SSO=1)) that is more stable method has no exponentially growing parts.

### How to use
```python
import get_TR from smatrix
```
Import and use `get_TR` from `smatrix.py`. 

It has 5 arguments (`layers`, `n_0`, `q_0`, `ws` and `sp`).

1. `layers`: Layer list. Unit of the thickness parameter in a layer is `cm`. 
   Please refer [example.ipynb](https://github.com/phykn/smatrix/blob/main/example.ipynb) for details.
2. `n_0`: Refractive index of the top layer.
3. `q_0`: Incidence angle of the top layer. (0 <= `q_0` < 90, unit:degrees)
4. `ws`: Wavernumber (unit: cm-1)
5. `sp`: Polarization (`'s'` or `'p'`) 
