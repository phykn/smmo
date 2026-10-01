import numpy as np
import pytest
import re
import subprocess
import sys

from pathlib import Path

from smmo import SMMO, make_config, make_layer


def constant_layer(n: float, thickness: float = 0.0, coherent: bool = False) -> dict:
    return layer(n, thickness=thickness, coherent=coherent)


def layer(
    n: float,
    k: float = 0.0,
    thickness: float = 0.0,
    coherent: bool = False,
) -> dict:
    return make_layer(
        n=np.array([n], dtype=float),
        k=np.array([k], dtype=float),
        thickness=thickness,
        coherent=coherent,
    )


def run_stack(layers: list[dict], incidence: float = 0.0, polarization: str = "s") -> dict:
    config = make_config(np.array([1000.0], dtype=float), incidence, polarization)
    return SMMO(layers, config)()


def test_air_glass_interface_returns_flux_transmittance() -> None:
    result = run_stack([
        constant_layer(1.0),
        constant_layer(1.5),
    ])

    assert set(result) == {"T", "R", "A"}
    np.testing.assert_allclose(result["T"], [0.96])
    np.testing.assert_allclose(result["R"], [0.04])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["T"] + result["R"], [1.0])


def test_lossless_air_glass_air_incoherent_slab_conserves_returned_intensity() -> None:
    result = run_stack([
        constant_layer(1.0),
        constant_layer(1.5, thickness=1000.0, coherent=False),
        constant_layer(1.0),
    ])

    np.testing.assert_allclose(result["T"], [12 / 13])
    np.testing.assert_allclose(result["R"], [1 / 13])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["T"] + result["R"], [1.0])


def test_brewster_angle_suppresses_p_polarized_reflection() -> None:
    brewster_angle = np.degrees(np.arctan(1.5 / 1.0))

    result = run_stack([
        constant_layer(1.0),
        constant_layer(1.5),
    ], incidence=brewster_angle, polarization="p")

    np.testing.assert_allclose(result["T"], [1.0])
    np.testing.assert_allclose(result["R"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)


def test_glass_air_interface_returns_flux_transmittance() -> None:
    result = run_stack([
        constant_layer(1.5),
        constant_layer(1.0),
    ])

    np.testing.assert_allclose(result["T"], [0.96])
    np.testing.assert_allclose(result["R"], [0.04])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["T"] + result["R"], [1.0])


def test_quarter_wave_coating_cancels_front_surface_reflection() -> None:
    n_coating = np.sqrt(1.5)
    thickness = (1.0 / 1000.0) / (4 * n_coating)

    result = run_stack([
        constant_layer(1.0),
        constant_layer(n_coating, thickness=thickness, coherent=True),
        constant_layer(1.5),
    ])

    np.testing.assert_allclose(result["R"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["T"], [1.0])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)


def test_matching_exit_medium_keeps_lossless_stack_intensity_sum_at_one() -> None:
    n_coating = np.sqrt(1.5)
    thickness = (1.0 / 1000.0) / (4 * n_coating)

    result = run_stack([
        constant_layer(1.0),
        constant_layer(n_coating, thickness=thickness, coherent=True),
        constant_layer(1.5),
        constant_layer(1.0),
    ])

    np.testing.assert_allclose(result["T"], [0.96])
    np.testing.assert_allclose(result["R"], [0.04])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["T"] + result["R"], [1.0])


def test_lossy_coherent_film_matches_single_film_solution() -> None:
    result = run_stack([
        constant_layer(1.0),
        layer(1.5, k=0.1, thickness=0.001, coherent=True),
        constant_layer(1.0),
    ])

    np.testing.assert_allclose(result["T"], [0.268621846195413])
    np.testing.assert_allclose(result["R"], [0.021741896060635])
    np.testing.assert_allclose(result["A"], [0.709636257743952])
    np.testing.assert_allclose(result["T"] + result["R"] + result["A"], [1.0])


def test_lossy_incoherent_slab_matches_intensity_sum_solution() -> None:
    result = run_stack([
        constant_layer(1.0),
        layer(1.5, k=0.1, thickness=0.001, coherent=False),
        constant_layer(1.0),
    ])

    np.testing.assert_allclose(result["T"], [0.262657558533408])
    np.testing.assert_allclose(result["R"], [0.0446383802595634])
    np.testing.assert_allclose(result["A"], [0.692704061207029])
    np.testing.assert_allclose(result["T"] + result["R"] + result["A"], [1.0])


def test_repeated_calls_are_stable_and_do_not_mutate_input_layers() -> None:
    layers = [
        constant_layer(1.0),
        constant_layer(2.0, thickness=0.1, coherent=True),
        constant_layer(1.5, thickness=10.0, coherent=False),
        constant_layer(1.0),
    ]
    original_length = len(layers)
    original = [
        {key: value.copy() if isinstance(value, np.ndarray) else value for key, value in layer.items()}
        for layer in layers
    ]
    config = make_config(np.array([1000.0], dtype=float), 0.0, "s")
    smmo = SMMO(layers, config)

    first = smmo()
    second = smmo()

    assert len(layers) == original_length
    assert len(smmo.layers) == original_length
    np.testing.assert_allclose(first["T"], second["T"])
    np.testing.assert_allclose(first["R"], second["R"])
    np.testing.assert_allclose(first["A"], second["A"])
    for before, after in zip(original, layers):
        assert before.keys() == after.keys()
        for key in before:
            np.testing.assert_array_equal(before[key], after[key])


def test_smatrix_component_calculation_does_not_mutate_input_layers() -> None:
    layers = [
        constant_layer(1.0, coherent=True),
        constant_layer(1.5, coherent=True),
    ]
    original_length = len(layers)
    config = make_config(np.array([1000.0], dtype=float), 0.0, "s")
    smmo = SMMO(layers, config)

    smmo.get_smatrix_components(layers)

    assert len(layers) == original_length


def test_invalid_polarization_is_rejected() -> None:
    config = make_config(np.array([1000.0], dtype=float), 0.0, "x")

    with pytest.raises(ValueError, match="polarization"):
        SMMO([constant_layer(1.0), constant_layer(1.5)], config)


def test_mismatched_spectral_lengths_are_rejected() -> None:
    config = make_config(np.array([1000.0, 1100.0], dtype=float), 0.0, "s")

    with pytest.raises(ValueError, match="size mismatch"):
        SMMO([constant_layer(1.0), constant_layer(1.5)], config)


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_total_internal_reflection_has_finite_power(polarization: str) -> None:
    with np.errstate(divide="raise", invalid="raise"):
        result = run_stack([
            constant_layer(1.5),
            constant_layer(1.0),
        ], incidence=60.0, polarization=polarization)

    np.testing.assert_allclose(result["T"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["R"], [1.0])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_absorbing_exit_interface_matches_fresnel_power(polarization: str) -> None:
    angle = np.deg2rad(50.0)
    n = 1.5 + 0.5j
    cos_i = np.cos(angle)
    cos_j = np.sqrt(1 - (np.sin(angle) / n) ** 2)
    if polarization == "s":
        r = (cos_i - n * cos_j) / (cos_i + n * cos_j)
        t = 2 * cos_i / (cos_i + n * cos_j)
        flux = np.real(n * cos_j) / cos_i
    else:
        r = (n * cos_i - cos_j) / (n * cos_i + cos_j)
        t = 2 * cos_i / (n * cos_i + cos_j)
        flux = np.real(n * np.conj(cos_j)) / cos_i

    result = run_stack([
        constant_layer(1.0),
        layer(1.5, k=0.5),
    ], incidence=50.0, polarization=polarization)

    np.testing.assert_allclose(result["T"], [flux * abs(t) ** 2])
    np.testing.assert_allclose(result["R"], [abs(r) ** 2])
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)


def test_invalid_config_is_rejected_in_optimized_python() -> None:
    script = """
import numpy as np
from smmo import SMMO, make_config, make_layer
layers = [make_layer(np.ones(1), np.zeros(1), 0.0, False)] * 2
for angle, pol in [(95.0, "s"), (0.0, "x")]:
    try:
        SMMO(layers, make_config(np.ones(1), angle, pol))
    except ValueError:
        pass
    else:
        raise RuntimeError("invalid configuration was accepted")
"""
    result = subprocess.run(
        [sys.executable, "-O", "-c", script],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("polarization", ["s", "p"])
@pytest.mark.parametrize("incidence", [np.degrees(np.arcsin(1 / 1.5)), 60.0])
@pytest.mark.parametrize("coherent", [False, True])
def test_evanescent_gap_matches_single_slab_solution(
    polarization: str, incidence: float, coherent: bool
) -> None:
    thickness = 0.0001
    with np.errstate(divide="raise", invalid="raise"):
        result = run_stack([
            constant_layer(1.5),
            constant_layer(1.0, thickness=thickness, coherent=coherent),
            constant_layer(1.5),
        ], incidence=incidence, polarization=polarization)

    if not coherent:
        transmission = 0.0
    else:
        cos_i = np.cos(np.deg2rad(incidence))
        cos_j = np.sqrt(1 - (1.5 * np.sin(np.deg2rad(incidence))) ** 2 + 0j)
        q_i = 1.5 * cos_i if polarization == "s" else cos_i / 1.5
        q_j = cos_j
        coef = 2 * np.pi * 1000.0 * thickness
        delta = coef * q_j
        b = -1j * coef if q_j == 0 else -1j * np.sin(delta) / q_j
        c = -1j * q_j * np.sin(delta)
        t = 2 * q_i / (2 * q_i * np.cos(delta) + q_i ** 2 * b + c)
        transmission = abs(t) ** 2
    np.testing.assert_allclose(result["T"], [transmission], atol=1e-12)
    np.testing.assert_allclose(result["R"], [1 - transmission], atol=1e-12)
    np.testing.assert_allclose(result["A"], [0.0], atol=1e-12)


@pytest.mark.parametrize("polarization", ["s", "p"])
@pytest.mark.parametrize("incidence", [0.0, 35.0, 60.0])
def test_coherent_spectrum_matches_characteristic_matrix(
    polarization: str, incidence: float
) -> None:
    w = np.array([600.0, 1000.0, 1800.0])
    n = np.array([1.5, 1.0, 2.0 + 0.05j, 1.3])
    thickness = [0.0, 0.00007, 0.0001, 0.0]
    layers = [
        make_layer(np.full_like(w, index.real), np.full_like(w, index.imag), d, True)
        for index, d in zip(n, thickness)
    ]
    angle = np.deg2rad(incidence)
    cos_q = np.sqrt(1 - (n[0] * np.sin(angle) / n) ** 2)
    q = n * cos_q if polarization == "s" else cos_q / n
    flux = np.real(n * (cos_q if polarization == "s" else np.conj(cos_q)))
    expected_t, expected_r = [], []
    for wn in w:
        matrix = np.eye(2, dtype=complex)
        for i in (1, 2):
            delta = 2 * np.pi * wn * n[i] * cos_q[i] * thickness[i]
            matrix = matrix @ np.array([
                [np.cos(delta), -1j * np.sin(delta) / q[i]],
                [-1j * q[i] * np.sin(delta), np.cos(delta)],
            ])
        a = matrix[0, 0] + matrix[0, 1] * q[-1]
        b = matrix[1, 0] + matrix[1, 1] * q[-1]
        r = (q[0] * a - b) / (q[0] * a + b)
        t = 2 * q[0] / (q[0] * a + b)
        if polarization == "p":
            t *= n[0] / n[-1]
        expected_t.append(flux[-1] / flux[0] * abs(t) ** 2)
        expected_r.append(abs(r) ** 2)

    result = SMMO(layers, make_config(w, incidence, polarization))()
    np.testing.assert_allclose(result["T"], expected_t, atol=1e-12)
    np.testing.assert_allclose(result["R"], expected_r, atol=1e-12)
    np.testing.assert_allclose(result["A"], 1 - np.array(expected_t) - expected_r, atol=1e-12)


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_opaque_film_stays_finite_and_matches_front_interface(polarization: str) -> None:
    n = 0.05 + 8j
    angle = np.deg2rad(60.0)
    cos_i = np.cos(angle)
    cos_j = np.sqrt(1 - (np.sin(angle) / n) ** 2)
    if polarization == "s":
        r = (cos_i - n * cos_j) / (cos_i + n * cos_j)
    else:
        r = (n * cos_i - cos_j) / (n * cos_i + cos_j)
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        result = run_stack([
            constant_layer(1.0),
            layer(n.real, n.imag, thickness=10.0, coherent=True),
            constant_layer(1.5),
        ], incidence=60.0, polarization=polarization)
    np.testing.assert_allclose(result["T"], [0.0], atol=1e-12)
    np.testing.assert_allclose(result["R"], [abs(r) ** 2])
    np.testing.assert_allclose(result["A"], [1 - abs(r) ** 2])


@pytest.mark.parametrize("count", [0, 1])
def test_missing_boundary_layers_are_rejected(count: int) -> None:
    with pytest.raises(ValueError, match="minimum number"):
        run_stack([constant_layer(1.0)] * count)


@pytest.mark.parametrize("incidence", [-1.0, 90.0, np.nan])
def test_invalid_incidence_is_rejected(incidence: float) -> None:
    with pytest.raises(ValueError, match="theta"):
        run_stack([constant_layer(1.0), constant_layer(1.5)], incidence=incidence)


@pytest.mark.parametrize("w", [np.array([]), np.ones((1, 1)), np.array(1000.0)])
def test_invalid_wavenumber_shape_is_rejected(w: np.ndarray) -> None:
    with pytest.raises(ValueError, match="one-dimensional"):
        SMMO([constant_layer(1.0), constant_layer(1.5)], make_config(w, 0.0, "s"))


@pytest.mark.parametrize("key", ["n", "k"])
def test_multidimensional_layer_data_is_rejected(key: str) -> None:
    layers = [constant_layer(1.0), constant_layer(1.5)]
    layers[1][key] = np.ones((1, 1))
    with pytest.raises(ValueError, match="same shape"):
        run_stack(layers)


def test_packaging_declares_supported_python_version() -> None:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = pyproject_path.read_text()
    match = re.search(r'^requires-python\s*=\s*"([^"]+)"$', pyproject, re.MULTILINE)

    assert match is not None
    assert match.group(1) == ">=3.9"
