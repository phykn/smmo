import numpy as np

from typing import Any


class SMMO:
    def __init__(
        self,
        layers: list[dict[str, Any]],
        config: dict[str, Any],
    ) -> None:
        self.check_data(layers, config)
        self.check_layers(layers)
        self.check_config(config)

        self.layers = layers
        self.n_0 = layers[0]["n"] + 1j * layers[0]["k"]
        self.wn = config["w"]
        self.theta0 = config["theta"]
        self.pol = config["pol"]

    @staticmethod
    def check_data(
        layers: list[dict[str, Any]],
        config: dict[str, Any],
    ) -> bool:
        lengths = [len(config["w"])]

        for layer in layers:
            lengths.append(len(layer["n"]))
            lengths.append(len(layer["k"]))

        assert len(np.unique(lengths)) == 1, "size mismatch"

        return True

    @staticmethod
    def check_layers(
        layers: list[dict[str, Any]]
    ) -> bool:
        assert len(layers) > 1, "The minimum number of layers is 2."

        return True

    @staticmethod
    def check_config(
        config: dict[str, Any]
    ) -> bool:
        assert 0 <= config["theta"] < 90, "theta must be set between 0 and 90."
        assert config["pol"] in ["s", "p"], 'polarization not in ["s", "p"]'

        return True

    @staticmethod
    def get_cos_qi(
        n_0: np.ndarray,
        n_i: np.ndarray,
        theta0: float,
    ) -> np.ndarray:
        theta0_arr = np.full(len(n_0), complex(theta0 * np.pi / 180, 0), dtype=complex)

        return np.sqrt(1 - np.square((n_0 * np.sin(theta0_arr) / n_i)))

    @staticmethod
    def get_kz(
        wavenumbers: np.ndarray,
        n_i: np.ndarray,
        cos_qi: np.ndarray,
    ) -> np.ndarray:
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
        if coeff == "r" and pol == "s":
            return (n_i * cos_qi - n_j * cos_qj) / (n_i * cos_qi + n_j * cos_qj)
        if coeff == "r" and pol == "p":
            return (n_j * cos_qi - n_i * cos_qj) / (n_j * cos_qi + n_i * cos_qj)
        if coeff == "t" and pol == "s":
            return 2 * n_i * cos_qi / (n_i * cos_qi + n_j * cos_qj)
        if coeff == "t" and pol == "p":
            return 2 * n_i * cos_qi / (n_j * cos_qi + n_i * cos_qj)

        raise ValueError("Invalid parameters for reflection/transmission coefficient calculation")

    def get_flux_factor(
        self,
        layers: list[dict[str, Any]],
    ) -> np.ndarray:
        n_i = layers[0]["n"] + 1j * layers[0]["k"]
        n_j = layers[-1]["n"] + 1j * layers[-1]["k"]
        cos_qi = self.get_cos_qi(self.n_0, n_i, self.theta0)
        cos_qj = self.get_cos_qi(self.n_0, n_j, self.theta0)

        return np.real(n_j * cos_qj) / np.real(n_i * cos_qi)

    def get_smatrix_components(
        self,
        layers: list[dict[str, Any]],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        layer_f = layers[-1].copy()
        layer_f["thickness"] = 0.0
        layer_f["coherent"] = True
        layers = [*layers, layer_f]

        num_layer = len(layers)
        num_data = len(self.wn)

        n = np.zeros((num_layer, num_data), dtype=complex)
        cos_q = np.zeros((num_layer, num_data), dtype=complex)
        phase = np.zeros((num_layer, num_data), dtype=complex)

        for i, layer in enumerate(layers):
            n[i] = layer["n"] + 1j * layer["k"]
            cos_q[i] = self.get_cos_qi(self.n_0, n[i], self.theta0)
            phase[i] = np.exp(
                1j * self.get_kz(self.wn, n[i], cos_q[i]) * layer["thickness"]
            )

        s_11 = np.full(num_data, 1, dtype=complex)
        s_12 = np.full(num_data, 0, dtype=complex)
        s_21 = np.full(num_data, 0, dtype=complex)
        s_22 = np.full(num_data, 1, dtype=complex)

        for i in range(num_layer - 1):
            j = i + 1
            t = self.get_fresnel_coeff_ij(n[i], n[j], cos_q[i], cos_q[j], coeff="t", pol=self.pol)
            r = self.get_fresnel_coeff_ij(n[i], n[j], cos_q[i], cos_q[j], coeff="r", pol=self.pol)

            i_11 = 1 / t
            i_12 = r / t
            i_21 = i_12
            i_22 = i_11

            s_11 = phase[i] * s_11 / (i_11 - phase[i] * s_12 * i_21)
            s_12 = ((phase[i] * s_12 * i_22 - i_12) * phase[j]) / (i_11 - phase[i] * s_12 * i_21)
            s_21 = s_22 * i_21 * s_11 + s_21
            s_22 = s_22 * i_21 * s_12 + s_22 * i_22 * phase[j]

        return s_11, s_12, s_21, s_22

    def get_tr_matrix_components(
        self,
        layers: list[dict[str, Any]],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        s_11, _, s_21, _ = self.get_smatrix_components(layers)
        t_12 = self.get_flux_factor(layers) * (s_11 * np.conjugate(s_11)).real
        r_12 = (s_21 * np.conjugate(s_21)).real

        s_11_r, _, s_21_r, _ = self.get_smatrix_components(layers[::-1])
        t_21 = self.get_flux_factor(layers[::-1]) * (s_11_r * np.conjugate(s_11_r)).real
        r_21 = (s_21_r * np.conjugate(s_21_r)).real

        return t_12, r_12, t_21, r_21

    @staticmethod
    def expand_layers(
        layers: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        pseudo_layers = []

        for layer in layers:
            if layer["coherent"]:
                pseudo_layers.append(layer)
            else:
                pseudo_layer = {
                    "n": layer["n"],
                    "k": layer["k"],
                    "thickness": 0.0,
                    "coherent": True,
                }
                pseudo_layers.append(pseudo_layer)
                pseudo_layers.append(layer)
                pseudo_layers.append(pseudo_layer)

        return pseudo_layers

    @staticmethod
    def split_blocks(
        layers: list[dict[str, Any]]
    ) -> list[list[dict[str, Any]]]:
        blocks = []
        block = layers[0:1]

        for layer in layers[1:]:
            if block[-1]["coherent"] == layer["coherent"]:
                block.append(layer)
            else:
                blocks.append(block)
                block = [layer]

        blocks.append(block)

        return blocks

    def __call__(
        self
    ) -> dict[str, np.ndarray]:
        pseudo_layers = self.expand_layers(self.layers)
        blocks = self.split_blocks(pseudo_layers)

        t_12_block = []
        r_12_block = []
        t_21_block = []
        r_21_block = []

        num_data = len(self.wn)
        t_12_total = np.full((2, num_data), 1, dtype=float)
        t_21_total = np.full((2, num_data), 1, dtype=float)
        r_12_total = np.full((2, num_data), 0, dtype=float)
        r_21_total = np.full((2, num_data), 0, dtype=float)

        for b in blocks:
            tr = self.get_tr_matrix_components(b)
            t_12_block.append(tr[0])
            r_12_block.append(tr[1])
            t_21_block.append(tr[2])
            r_21_block.append(tr[3])

        t_12_arr = np.array(t_12_block)
        r_12_arr = np.array(r_12_block)
        t_21_arr = np.array(t_21_block)
        r_21_arr = np.array(r_21_block)

        for i in range(len(blocks)):
            norm = 1 - r_21_total[0] * r_12_arr[i]
            t_12_total[1] = t_12_total[0] * t_12_arr[i] / norm
            r_12_total[1] = r_12_total[0] + t_12_total[0] * t_21_total[0] * r_12_arr[i] / norm
            t_21_total[1] = t_21_total[0] * t_21_arr[i] / norm
            r_21_total[1] = r_21_arr[i] + t_12_arr[i] * t_21_arr[i] * r_21_total[0] / norm

            t_12_total[0] = t_12_total[1]
            r_12_total[0] = r_12_total[1]
            t_21_total[0] = t_21_total[1]
            r_21_total[0] = r_21_total[1]

        absorption = 1 - t_12_total[1] - r_12_total[1]
        absorption = np.where(np.isclose(absorption, 0.0, atol=1e-12), 0.0, absorption)

        return {
            "T": t_12_total[1],
            "R": r_12_total[1],
            "A": absorption,
        }
