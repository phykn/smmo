import numpy as np

from .data import Config, Layer


class SMMO:
    """Calculate spectral power for a stack ordered from entrance to exit."""

    def __init__(self, layers: list[Layer], config: Config) -> None:
        self.check_layers(layers)
        self.check_config(config)
        self.check_data(layers, config)

        self.layers = layers
        self.n_0 = layers[0]["n"] + 1j * layers[0]["k"]
        self.wn = config["w"]
        self.theta0 = config["theta"]
        self.pol = config["pol"]

    def __call__(self) -> dict[str, np.ndarray]:
        blocks = self.split_blocks(self.expand_layers(self.layers))
        t_12 = np.ones(len(self.wn))
        t_21 = np.ones(len(self.wn))
        r_12 = np.zeros(len(self.wn))
        r_21 = np.zeros(len(self.wn))

        for block in blocks:
            t_12, r_12, t_21, r_21 = self._cascade(
                (t_12, r_12, t_21, r_21), self.get_tr_matrix_components(block)
            )

        absorption = 1 - t_12 - r_12
        absorption = np.where(np.isclose(absorption, 0.0, atol=1e-12), 0.0, absorption)
        return {"T": t_12, "R": r_12, "A": absorption}

    @staticmethod
    def _cascade(
        left: tuple[np.ndarray, ...], right: tuple[np.ndarray, ...]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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

    @staticmethod
    def check_layers(layers: list[Layer]) -> bool:
        if len(layers) < 2:
            raise ValueError("The minimum number of layers is 2.")
        return True

    @staticmethod
    def check_config(config: Config) -> bool:
        if not 0 <= config["theta"] < 90:
            raise ValueError("theta must be set between 0 and 90.")
        if config["pol"] not in ("s", "p"):
            raise ValueError('polarization not in ["s", "p"]')
        return True

    @staticmethod
    def check_data(layers: list[Layer], config: Config) -> bool:
        w = config["w"]
        if np.ndim(w) != 1 or len(w) == 0:
            raise ValueError("wavenumber must be a nonempty one-dimensional array")
        for layer in layers:
            for key in ("n", "k"):
                if np.ndim(layer[key]) != 1 or len(layer[key]) != len(w):
                    raise ValueError("size mismatch: n and k must have the same shape as wavenumber")
        return True

    @staticmethod
    def expand_layers(layers: list[Layer]) -> list[Layer]:
        """Surround each incoherent layer with zero-thickness coherent boundaries."""
        expanded = []
        for layer in layers:
            if layer["coherent"]:
                expanded.append(layer)
            else:
                boundary: Layer = {
                    "n": layer["n"],
                    "k": layer["k"],
                    "thickness": 0.0,
                    "coherent": True,
                }
                expanded.extend([boundary, layer, boundary])
        return expanded

    @staticmethod
    def split_blocks(layers: list[Layer]) -> list[list[Layer]]:
        blocks = []
        block = layers[:1]
        for layer in layers[1:]:
            if block[-1]["coherent"] == layer["coherent"]:
                block.append(layer)
            else:
                blocks.append(block)
                block = [layer]
        blocks.append(block)
        return blocks

    def get_tr_matrix_components(
        self, layers: list[Layer]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        s_11, _, s_21, _ = self.get_smatrix_components(layers)
        t_12 = self.get_flux_factor(layers) * np.abs(s_11) ** 2
        r_12 = np.abs(s_21) ** 2

        reverse = layers[::-1]
        s_11, _, s_21, _ = self.get_smatrix_components(reverse)
        t_21 = self.get_flux_factor(reverse) * np.abs(s_11) ** 2
        r_21 = np.abs(s_21) ** 2
        return t_12, r_12, t_21, r_21

    def get_flux_factor(self, layers: list[Layer]) -> np.ndarray:
        n_i = layers[0]["n"] + 1j * layers[0]["k"]
        n_j = layers[-1]["n"] + 1j * layers[-1]["k"]
        cos_i = self.get_cos_qi(self.n_0, n_i, self.theta0)
        cos_j = self.get_cos_qi(self.n_0, n_j, self.theta0)
        if self.pol == "p":
            cos_i, cos_j = np.conj(cos_i), np.conj(cos_j)
        flux_i = np.real(n_i * cos_i)
        flux_j = np.real(n_j * cos_j)

        # An evanescent port carries no normal power; equal zero-flux ports
        # retain the unit factor needed for propagation within one medium.
        factor = np.where(flux_j == 0, 1.0, 0.0)
        return np.divide(flux_j, flux_i, out=factor, where=flux_i != 0)

    def get_smatrix_components(
        self, layers: list[Layer]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        n = np.array([layer["n"] + 1j * layer["k"] for layer in layers])
        cos_q = self.get_cos_qi(self.n_0, n, self.theta0)
        thickness = np.array([layer["thickness"] for layer in layers])
        phase = np.exp(1j * self.get_kz(self.wn, n, cos_q) * thickness[:, None])
        zero = np.zeros(len(self.wn), dtype=complex)
        if len(layers) == 1:
            return phase[0], zero.copy(), zero, np.ones(len(self.wn), dtype=complex)

        q = n * cos_q if self.pol == "s" else cos_q / n
        cos_0 = self.get_cos_qi(self.n_0, self.n_0, self.theta0)
        q_0 = self.n_0 * cos_0 if self.pol == "s" else cos_0 / self.n_0
        ref = np.where(q[0] != 0, q[0], q_0)
        total = self._interface(q[0], ref)

        # Embed each propagation step in the same nonzero-admittance medium.
        # This avoids singular interfaces at an internal layer's critical angle.
        for i in range(1, len(layers)):
            coef = 2 * np.pi * self.wn * thickness[i]
            if self.pol == "p":
                coef = coef * n[i] ** 2
            delta = coef * q[i]
            change = np.expm1(2j * delta)
            a = 1 + change / 2
            b = np.divide(-change / 2, q[i], out=-1j * coef.astype(complex), where=q[i] != 0)
            c = -change * q[i] / 2
            norm = 2 * ref * a + ref ** 2 * b + c
            t = 2 * ref * phase[i] / norm
            r = (ref ** 2 * b - c) / norm
            total = self._cascade(total, (t, r, t, r))

        t_f, r_f, t_b, r_b = self._cascade(total, self._interface(ref, q[-1]))
        uniform = (q[0] == 0) & np.all(n == n[0], axis=0)
        propagation = np.prod(phase[1:], axis=0)
        t_f = np.where(uniform, propagation, t_f)
        t_b = np.where(uniform, propagation, t_b)
        r_f = np.where(uniform, 0.0, r_f)
        r_b = np.where(uniform, 0.0, r_b)
        if self.pol == "p":
            t_f = t_f * n[0] / n[-1]
            t_b = t_b * n[-1] / n[0]

        # Preserve the entrance-phase convention of the amplitude components.
        return t_f * phase[0], r_b, r_f * phase[0], t_b

    @staticmethod
    def _interface(
        q_i: np.ndarray, q_j: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        norm = q_i + q_j
        r = (q_i - q_j) / norm
        return 2 * q_i / norm, r, 2 * q_j / norm, -r

    @staticmethod
    def get_cos_qi(n_0: np.ndarray, n_i: np.ndarray, theta0: float) -> np.ndarray:
        return np.sqrt(1 - (n_0 * np.sin(np.deg2rad(theta0)) / n_i) ** 2 + 0j)

    @staticmethod
    def get_kz(wavenumbers: np.ndarray, n_i: np.ndarray, cos_qi: np.ndarray) -> np.ndarray:
        return 2 * np.pi * n_i * wavenumbers * cos_qi

    @staticmethod
    def get_fresnel_coeff_ij(
        n_i: np.ndarray,
        n_j: np.ndarray,
        cos_qi: np.ndarray,
        cos_qj: np.ndarray,
        coeff: str = "r",
        pol: str = "s",
    ) -> np.ndarray:
        if pol == "s":
            a, b = n_i * cos_qi, n_j * cos_qj
        elif pol == "p":
            a, b = n_j * cos_qi, n_i * cos_qj
        else:
            raise ValueError("Invalid parameters for reflection/transmission coefficient calculation")

        if coeff == "r":
            return (a - b) / (a + b)
        if coeff == "t":
            return 2 * n_i * cos_qi / (a + b)
        raise ValueError("Invalid parameters for reflection/transmission coefficient calculation")
