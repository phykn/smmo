import numpy as np

from typing import Any, Dict


def make_config(
    wavenumber: np.ndarray,
    incidence: float,
    polarization: str
) -> Dict[str, Any]:
    return dict(
        w = wavenumber,
        q = incidence,
        p = polarization
    )


def make_layer(
    n: np.ndarray,
    k: np.ndarray,
    thickness: float,
    coherence: bool
) -> Dict[str, Any]:
    return dict(
        n = n,
        k = k,
        thickness = thickness,
        coherence = coherence
    )