### Scattering-Matrix method for multilayer optics
Transfer-matrix method(TMM) is a simple, accurate and fast method used to understand the propagation of electromagnetic waves in a multilayer structure. However, TMM is numerically unstable when a thick layer in the structure has high absorption because exponential parts cause overflow. To overcome this problem, scattering-matrix method(SMM) was presented (PRB-1988, JOSAA-1988, SPIE-2010) which is more stable method has no exponentially growing parts.

### How to use
Import get_TR from smatrix.py. It has 4 arguments (layers, n_0, q_0, ws and sp)

1. layers: layers is list of layers. Please refer the example.ipynb
2. n_0: Refractive index of the top layer.
3. q_0: Incidence angle of the top layer.
4. ws: Wavernumber (cm-1)
5. sp: Polarization ('s' or 'p') 
