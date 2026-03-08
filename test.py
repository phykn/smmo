import numpy as np

from smmo import SMMO, make_config, make_layer


W_VALS = np.linspace(500, 2000, 100)


def run_test(
    name: str,
    layers: list[dict],
    cfg: dict,
) -> None:
    res = SMMO(layers, cfg)()
    t_avg = np.mean(res["T"])
    r_avg = np.mean(res["R"])
    print(f"[{name}] T_avg: {t_avg:.4f}, R_avg: {r_avg:.4f}, Sum: {t_avg + r_avg:.4f}")


def test_vacuum() -> None:
    cfg = make_config(W_VALS, 0, "s")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
    ]
    run_test("Vacuum", layers, cfg)


def test_air_glass() -> None:
    cfg = make_config(W_VALS, 0, "s")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 1.5), np.zeros(100), 0, False),
    ]
    run_test("Air-Glass Interface", layers, cfg)


def test_thick_slab() -> None:
    cfg = make_config(W_VALS, 0, "s")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 1.5), np.zeros(100), 1000.0, False),
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
    ]
    run_test("Thick Incoherent Slab", layers, cfg)


def test_thin_film() -> None:
    cfg = make_config(W_VALS, 0, "s")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 1.5), np.zeros(100), 0.5, True),
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
    ]
    run_test("Thin Coherent Film", layers, cfg)


def test_lossy_material() -> None:
    cfg = make_config(W_VALS, 0, "s")
    k_vals = np.full(100, 0.1)
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 1.5), k_vals, 1.0, True),
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
    ]
    run_test("Lossy Material", layers, cfg)


def test_brewster_angle() -> None:
    # Brewster angle for n1=1, n2=1.5 is arctan(1.5/1) ~ 56.31 degrees
    cfg = make_config(W_VALS, 56.31, "p")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 1.5), np.zeros(100), 0, False),
    ]
    run_test("Brewster Angle (p-pol)", layers, cfg)


def test_anti_reflection() -> None:
    # n_layer = sqrt(1 * 1.5) ~ 1.225
    # thickness = lambda / (4 * n)
    w_mid = 1000.0
    n_ar = np.sqrt(1.5)
    d_ar = (1.0 / w_mid) / (4 * n_ar)
    cfg = make_config(W_VALS, 0, "s")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, n_ar), np.zeros(100), d_ar, True),
        make_layer(np.full(100, 1.5), np.zeros(100), 0, False),
    ]
    run_test("Anti-Reflection Coating", layers, cfg)


def test_bragg_mirror() -> None:
    n_l, n_h = 1.38, 2.35
    w_mid = 1000.0
    d_l = (1.0 / w_mid) / (4 * n_l)
    d_h = (1.0 / w_mid) / (4 * n_h)
    
    cfg = make_config(W_VALS, 0, "s")
    layers = [make_layer(np.full(100, 1.0), np.zeros(100), 0, False)]
    for _ in range(5):
        layers.append(make_layer(np.full(100, n_h), np.zeros(100), d_h, True))
        layers.append(make_layer(np.full(100, n_l), np.zeros(100), d_l, True))
    layers.append(make_layer(np.full(100, 1.5), np.zeros(100), 0, False))
    
    run_test("Bragg Mirror (5 pairs)", layers, cfg)


def test_mixed_coherence() -> None:
    cfg = make_config(W_VALS, 0, "s")
    layers = [
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
        make_layer(np.full(100, 2.0), np.zeros(100), 0.1, True),
        make_layer(np.full(100, 1.5), np.zeros(100), 10.0, False),
        make_layer(np.full(100, 1.0), np.zeros(100), 0, False),
    ]
    run_test("Mixed Coherence", layers, cfg)


def test_wide_spectrum() -> None:
    w_wide = np.linspace(100, 5000, 500)
    cfg = make_config(w_wide, 30, "s")
    layers = [
        make_layer(np.full(500, 1.0), np.zeros(500), 0, False),
        make_layer(np.full(500, 1.5), np.zeros(500), 1.0, True),
        make_layer(np.full(500, 1.0), np.zeros(500), 0, False),
    ]
    run_test("Wide Spectrum Sweep", layers, cfg)


def main() -> None:
    test_vacuum()
    test_air_glass()
    test_thick_slab()
    test_thin_film()
    test_lossy_material()
    test_brewster_angle()
    test_anti_reflection()
    test_bragg_mirror()
    test_mixed_coherence()
    test_wide_spectrum()


if __name__ == "__main__":
    main()
