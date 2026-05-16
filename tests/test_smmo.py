import numpy as np
import pytest
import re

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
    config = make_config(np.array([1000.0], dtype=float), 0.0, "s")
    smmo = SMMO(layers, config)

    first = smmo()
    second = smmo()

    assert len(layers) == original_length
    assert len(smmo.layers) == original_length
    np.testing.assert_allclose(first["T"], second["T"])
    np.testing.assert_allclose(first["R"], second["R"])
    np.testing.assert_allclose(first["A"], second["A"])


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

    with pytest.raises(AssertionError, match="polarization"):
        SMMO([constant_layer(1.0), constant_layer(1.5)], config)


def test_mismatched_spectral_lengths_are_rejected() -> None:
    config = make_config(np.array([1000.0, 1100.0], dtype=float), 0.0, "s")

    with pytest.raises(AssertionError, match="size mismatch"):
        SMMO([constant_layer(1.0), constant_layer(1.5)], config)


def test_packaging_declares_supported_python_version() -> None:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = pyproject_path.read_text()
    match = re.search(r'^requires-python\s*=\s*"([^"]+)"$', pyproject, re.MULTILINE)

    assert match is not None
    assert match.group(1) == ">=3.9"
