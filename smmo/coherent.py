"""Scattering composition and bidirectional power for a coherent optical block."""

import numpy as np

from .data import Config, Layer


Coefficients = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def cascade(left: Coefficients, right: Coefficients) -> Coefficients:
    """Join amplitude or power coefficients ordered as T/R forward, T/R backward."""
    t_f, r_f, t_b, r_b = left
    u_f, v_f, u_b, v_b = right
    norm = 1 - r_b * v_f
    # Perfect reflectors disconnect the path between the two blocks.
    coupling = np.divide(1, norm, out=np.zeros_like(norm), where=norm != 0)
    return (
        t_f * u_f * coupling,
        r_f + t_f * t_b * v_f * coupling,
        t_b * u_b * coupling,
        v_b + u_f * u_b * r_b * coupling,
    )


def power_coefficients(layers: list[Layer], config: Config, n_0: np.ndarray) -> Coefficients:
    """Solve both directions with the same conserved transverse wavevector."""
    w, theta, pol = config["w"], config["theta"], config["pol"]
    n = np.array([layer["n"] + 1j * layer["k"] for layer in layers])
    thickness = np.array([layer["thickness"] for layer in layers])
    cos_q = np.sqrt(1 - (n_0 * np.sin(np.deg2rad(theta)) / n) ** 2 + 0j)
    phase = np.exp(2j * np.pi * w * n * cos_q * thickness[:, None])
    if len(layers) == 1:
        transmission = np.abs(phase[0]) ** 2
        zero = np.zeros(len(w))
        return transmission, zero, transmission, zero

    q = n * cos_q if pol == "s" else cos_q / n
    cos_0 = np.sqrt(1 - np.sin(np.deg2rad(theta)) ** 2 + 0j)
    q_0 = n_0 * cos_0 if pol == "s" else cos_0 / n_0
    ref = np.where(q[0] != 0, q[0], q_0)
    total = _interface(q[0], ref)

    # Embed internal propagation in a nonzero-admittance medium so the
    # critical-angle limit does not require a singular Fresnel interface.
    for i in range(1, len(layers) - 1):
        coef = 2 * np.pi * w * thickness[i]
        if pol == "p":
            coef = coef * n[i] ** 2
        delta = coef * q[i]
        change = np.expm1(2j * delta)
        a = 1 + change / 2
        b = np.divide(-change / 2, q[i], out=-1j * coef.astype(complex), where=q[i] != 0)
        c = -change * q[i] / 2
        norm = 2 * ref * a + ref ** 2 * b + c
        t = 2 * ref * phase[i] / norm
        r = (ref ** 2 * b - c) / norm
        total = cascade(total, (t, r, t, r))

    t_f, r_f, t_b, r_b = cascade(total, _interface(ref, q[-1]))
    uniform = (q[0] == 0) & np.all(n == n[0], axis=0)
    propagation = np.prod(phase[1:-1], axis=0)
    t_f = np.where(uniform, propagation, t_f)
    t_b = np.where(uniform, propagation, t_b)
    r_f = np.where(uniform, 0.0, r_f)
    r_b = np.where(uniform, 0.0, r_b)
    if pol == "p":
        t_f = t_f * n[0] / n[-1]
        t_b = t_b * n[-1] / n[0]

    # Preserve the boundary propagation convention for finite entrance/exit media.
    t_f = t_f * phase[0] * phase[-1]
    t_b = t_b * phase[0] * phase[-1]
    r_f = r_f * phase[0]
    r_b = r_b * phase[-1]
    flux = np.real(n * (cos_q if pol == "s" else np.conj(cos_q)))
    return (
        _flux_ratio(flux[0], flux[-1]) * np.abs(t_f) ** 2,
        np.abs(r_f) ** 2,
        _flux_ratio(flux[-1], flux[0]) * np.abs(t_b) ** 2,
        np.abs(r_b) ** 2,
    )


def _interface(q_i: np.ndarray, q_j: np.ndarray) -> Coefficients:
    norm = q_i + q_j
    r = (q_i - q_j) / norm
    return 2 * q_i / norm, r, 2 * q_j / norm, -r


def _flux_ratio(flux_i: np.ndarray, flux_j: np.ndarray) -> np.ndarray:
    # An evanescent port carries no normal power. Equal zero-flux ports retain
    # the unit factor needed for propagation within one medium.
    factor = np.where(flux_j == 0, 1.0, 0.0)
    return np.divide(flux_j, flux_i, out=factor, where=flux_i != 0)
