from collections.abc import Iterator

import numpy as np

from .coherent import cascade, power_coefficients
from .data import Config, Layer, validate


class SMMO:
    """Calculate spectral power using the current entrance-to-exit stack and config."""

    def __init__(self, layers: list[Layer], config: Config) -> None:
        validate(layers, config)
        self.layers = layers
        self.config = config

    def __call__(self) -> dict[str, np.ndarray]:
        validate(self.layers, self.config)
        n_0 = self.layers[0]["n"] + 1j * self.layers[0]["k"]
        one = np.ones(len(self.config["w"]))
        zero = np.zeros(len(self.config["w"]))
        total = one, zero, one, zero
        for block in self._blocks():
            total = cascade(total, power_coefficients(block, self.config, n_0))

        t, r, _, _ = total
        absorption = 1 - t - r
        absorption = np.where(np.isclose(absorption, 0.0, atol=1e-12), 0.0, absorption)
        return {"T": t, "R": r, "A": absorption}

    def _blocks(self) -> Iterator[list[Layer]]:
        """Yield coherent runs and isolated incoherent propagation layers."""
        block = []
        for layer in self.layers:
            if layer["coherent"]:
                block.append(layer)
            else:
                boundary: Layer = {**layer, "thickness": 0.0, "coherent": True}
                block.append(boundary)
                yield block
                yield [layer]
                block = [boundary]
        if block:
            yield block
