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


def validate(layers: list[Layer], config: Config) -> None:
    """Check the current stack before evaluating any optical expressions."""
    if len(layers) < 2:
        raise ValueError("The minimum number of layers is 2.")
    if not 0 <= config["theta"] < 90:
        raise ValueError("theta must be set between 0 and 90.")
    if config["pol"] not in ("s", "p"):
        raise ValueError('polarization not in ["s", "p"]')

    w = config["w"]
    if np.ndim(w) != 1 or len(w) == 0:
        raise ValueError("wavenumber must be a nonempty one-dimensional array")
    if not np.isrealobj(w) or not np.all(np.isfinite(w)):
        raise ValueError("wavenumber must contain finite real values")
    for layer in layers:
        for key in ("n", "k"):
            values = layer[key]
            if np.ndim(values) != 1 or len(values) != len(w):
                raise ValueError("size mismatch: n and k must have the same shape as wavenumber")
            if not np.isrealobj(values) or not np.all(np.isfinite(values)):
                raise ValueError(f"{key} must contain finite real values")
        d = layer["thickness"]
        if np.ndim(d) != 0 or not np.isrealobj(d) or not np.isfinite(d) or d < 0:
            raise ValueError("thickness must be a finite nonnegative real scalar")
