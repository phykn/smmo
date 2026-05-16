# SMMO

SMMO is a small Python package for multilayer optics. It calculates how much light is transmitted, reflected, and absorbed when light enters a stack of materials such as air, glass, coatings, films, and substrates.

The calculation uses the scattering-matrix method. This is useful for multilayer structures because it is more numerically stable than a plain transfer-matrix multiplication when layers are thick, lossy, or mixed between coherent and incoherent behavior.

## Purpose

SMMO is meant for one-dimensional optical stack calculations. It answers questions such as:

- How much light passes through this multilayer structure?
- How much light is reflected at the entrance side?
- How much light is absorbed inside lossy materials?
- How does the answer change across a spectrum, angle, polarization, or layer design?

The typical input is a list of layers with refractive index `n`, extinction coefficient `k`, thickness, and coherence behavior. The typical output is spectral `T`, `R`, and `A`, which stand for transmittance, reflectance, and absorptance.

## Why Use SMMO

Use SMMO when you need a direct numerical model of a multilayer optical structure instead of a single-interface Fresnel calculation. A single air-glass boundary can be calculated by hand, but real samples often include films, coatings, substrates, absorbing materials, or mixed coherent and incoherent regions.

SMMO keeps those pieces in one layer list. That makes it useful for checking whether a proposed stack should transmit, reflect, or absorb light before comparing against measured spectra or changing a design.

## When To Use It

SMMO is appropriate for planar multilayer structures where each layer is laterally uniform and can be described by `n`, `k`, and thickness. Common use cases include thin-film coatings, semiconductor films, optical substrates, anti-reflection coatings, Bragg mirrors, and spectrum fitting workflows.

It is also useful when the stack contains both thin coherent layers and thick incoherent layers. For example, a thin coating on a thick substrate can be modeled without treating the whole sample as fully coherent.

SMMO is not a full 3D electromagnetic solver. It does not model lateral patterns, surface roughness in detail, scattering from particles, diffraction gratings, or finite-size beam effects.

## Advantages

- **Stable multilayer calculation**: Uses scattering matrices to reduce numerical problems that can appear in direct transfer-matrix multiplication.
- **Mixed coherence support**: Allows coherent thin films and incoherent bulk layers in the same stack.
- **Loss handling**: Supports materials with `k > 0` and reports absorption through `A = 1 - T - R`.
- **Spectral calculation**: Accepts NumPy arrays, so one call can evaluate many wavenumber points.
- **Physical power transmittance**: Reports flux-corrected `T`, not just the squared transmission amplitude.
- **Small API**: The main workflow is `make_config(...)`, `make_layer(...)`, and `SMMO(layers, cfg)()`.

## Installation

```bash
pip install py-smmo
```

The PyPI distribution name is `py-smmo`. The Python import name remains `smmo`.

The package depends on NumPy. That dependency is declared in `pyproject.toml`.

## Quick Start

This example sends normally incident s-polarized light from air into a thin glass film and then back into air.

```python
import numpy as np

from smmo import SMMO, make_config, make_layer

w = np.linspace(500.0, 2000.0, 100)  # wavenumber in cm^-1

cfg = make_config(
    wavenumber=w,
    incidence=0.0,
    polarization="s",
)

layers = [
    make_layer(
        n=np.full_like(w, 1.0),
        k=np.zeros_like(w),
        thickness=0.0,
        coherent=False,
    ),
    make_layer(
        n=np.full_like(w, 1.5),
        k=np.zeros_like(w),
        thickness=0.01,
        coherent=True,
    ),
    make_layer(
        n=np.full_like(w, 1.0),
        k=np.zeros_like(w),
        thickness=0.0,
        coherent=False,
    ),
]

result = SMMO(layers, cfg)()

print(result["T"].mean())
print(result["R"].mean())
print(result["A"].mean())
```

`result["T"]`, `result["R"]`, and `result["A"]` are NumPy arrays with the same length as `w`.

## What The Inputs Mean

### Light Configuration

`make_config(...)` describes the incoming light.

| Argument | Meaning | Unit or value |
| --- | --- | --- |
| `wavenumber` | Spectral positions where the stack is evaluated. A larger array gives a spectrum. | `cm^-1` |
| `incidence` | Angle of the incoming light measured from the surface normal. `0` means straight into the stack. | degrees |
| `polarization` | Light polarization. Use `"s"` or `"p"`. | string |

The wavenumber array controls the length of every output array. For example, if `wavenumber` has 100 values, then `T`, `R`, and `A` each have 100 values.

### Layers

`make_layer(...)` describes one material layer. The order of the list matters: the first layer is where light comes from, and the last layer is where transmitted light exits.

| Argument | Meaning | Unit or value |
| --- | --- | --- |
| `n` | Real refractive index. This controls phase velocity and reflection at interfaces. | array |
| `k` | Extinction coefficient. This is the imaginary part of the refractive index and represents loss. | array |
| `thickness` | Physical layer thickness. Boundary media such as air can use `0.0`. | `cm` |
| `coherent` | Whether phase interference inside this layer is kept. | boolean |

The optical index used by the code is `N = n + i k`. If `k` is zero, the material itself has no absorption. If `k` is positive, some light can be absorbed inside the layer.

`n` and `k` should be arrays with the same length as `wavenumber`. Constant materials can be written with `np.full_like(w, value)` and `np.zeros_like(w)`.

## Coherent vs Incoherent Layers

Use `coherent=True` when phase matters. This is usually the right choice for thin films, coatings, and layers where interference fringes are physically meaningful.

Use `coherent=False` when the layer is thick enough, rough enough, or experimentally averaged enough that the phase is effectively washed out. This is usually used for bulk substrates or entrance and exit media.

Changing `coherent` changes the physical approximation. It is not just a performance option. A thin film marked as incoherent can give a different result because interference is intentionally removed.

## What The Outputs Mean

| Output | Meaning |
| --- | --- |
| `T` | Flux transmittance. This is the fraction of incoming optical power that exits the final layer. |
| `R` | Reflectance. This is the fraction of incoming optical power reflected back into the first layer. |
| `A` | Absorptance. This is computed as `1 - T - R`. |

For lossless stacks, `T + R = 1` and `A = 0`.

For lossy stacks, `T + R < 1` and the missing part is absorption, so `T + R + A = 1`.

`T` is not just `abs(t)^2`. When the entrance and exit media are different, the transmitted amplitude must be corrected by the optical flux factor. SMMO applies that correction internally, so `T` is the physical power transmittance.

## Examples

### Air To Glass Interface

At normal incidence, an air-to-glass boundary with `n = 1.5` reflects about 4% of incoming power and transmits about 96%. There is no absorption because `k = 0`.

```python
import numpy as np

from smmo import SMMO, make_config, make_layer

w = np.array([1000.0])
cfg = make_config(w, 0.0, "s")

layers = [
    make_layer(np.array([1.0]), np.array([0.0]), 0.0, False),
    make_layer(np.array([1.5]), np.array([0.0]), 0.0, False),
]

result = SMMO(layers, cfg)()

print(result["T"][0])  # about 0.96
print(result["R"][0])  # about 0.04
print(result["A"][0])  # about 0.00
```

This is a useful sanity check: because the stack is lossless, `T + R` should be `1`.

### Lossy Film

If a layer has `k > 0`, some light is absorbed. In that case, `T + R` is less than `1`, and `A` reports the absorbed fraction.

```python
import numpy as np

from smmo import SMMO, make_config, make_layer

w = np.array([1000.0])
cfg = make_config(w, 0.0, "s")

layers = [
    make_layer(np.array([1.0]), np.array([0.0]), 0.0, False),
    make_layer(np.array([1.5]), np.array([0.1]), 0.001, True),
    make_layer(np.array([1.0]), np.array([0.0]), 0.0, False),
]

result = SMMO(layers, cfg)()

print(result["T"][0])
print(result["R"][0])
print(result["A"][0])
print(result["T"][0] + result["R"][0] + result["A"][0])  # 1.0
```

Here `A` is positive because the middle layer has loss.

## Common Mistakes

Use centimeters for thickness. The code expects `thickness` in `cm`, so `1 mm` should be written as `0.1`, and `10 um` should be written as `0.001`.

Do not confuse `k` with an absorption coefficient. `k` is the extinction coefficient in the complex refractive index `n + i k`.

Make sure all spectral arrays have the same length. `wavenumber`, each layer's `n`, and each layer's `k` should line up point by point.

Remember that layer order is physical. `[air, film, glass]` and `[glass, film, air]` are different optical problems.

## Tests

Run the test suite with:

```bash
python -m pytest -q
```

Run the executable examples with:

```bash
python tests/test_examples.py
```

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
