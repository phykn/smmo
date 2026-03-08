import numpy as np

from typing import Any


def make_config(
    wavenumber: np.ndarray,
    incidence: float,
    polarization: str,
) -> dict[str, Any]:
    return {
        "w": wavenumber,
        "q": incidence,
        "p": polarization,
    }


def make_layer(
    n: np.ndarray,
    k: np.ndarray,
    thickness: float,
    coherence: bool,
) -> dict[str, Any]:
    return {
        "n": n,
        "k": k,
        "thickness": thickness,
        "coherence": coherence,
    }