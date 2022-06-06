import numpy as np
from typing import Dict, Any


def make_config(
    wavenumber: np.ndarray,
    incidence: float,
    polarization: str
) -> Dict[str, Any]:
    '''
    config module
    wavenumber: cm-1 (unit)
    incidence: degrees (unit)
    polarization: "s" or "p"
    '''
    return dict(
        w=wavenumber,
        q=incidence,
        p=polarization
    )


def make_layer(
    n: np.ndarray,
    k: np.ndarray,
    thickness: float,
    coherence: bool
) -> Dict[str, Any]:
    '''
    layer module
    n: refractive index
    k: absorption coeffeicient
    thickness: layer thickness, cm (unit)
    coherence: coherence (True) or incoherence (False)
    '''
    return dict(
        n=n,
        k=k,
        thickness=thickness,
        coherence=coherence
    )