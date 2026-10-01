from typing import TypedDict

import numpy as np


class Config(TypedDict):
    w: np.ndarray
    theta: float
    pol: str


class Layer(TypedDict):
    n: np.ndarray
    k: np.ndarray
    thickness: float
    coherent: bool


def make_config(
    wavenumber: np.ndarray,
    incidence: float,
    polarization: str,
) -> Config:
    """Describe light with wavenumbers in cm^-1 and incidence in degrees."""
    return {
        "w": wavenumber,
        "theta": incidence,
        "pol": polarization,
    }


def make_layer(
    n: np.ndarray,
    k: np.ndarray,
    thickness: float,
    coherent: bool,
) -> Layer:
    """Describe a material with index n + i k and thickness in cm."""
    return {
        "n": n,
        "k": k,
        "thickness": thickness,
        "coherent": coherent,
    }
