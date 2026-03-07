import numpy as np

from smmo import SMMO, make_config, make_layer


def main() -> None:
    config = make_config(
        wavenumber=np.arange(0, 10000, step=1000),
        incidence=0,
        polarization="s"
    )

    layers = [
        make_layer(
            n=np.full(10, 1),
            k=np.full(10, 0),
            thickness=0,
            coherence=False
        ),
        make_layer(
            n=np.full(10, 1.5),
            k=np.full(10, 0),
            thickness=0.01,
            coherence=True
        ),
        make_layer(
            n=np.full(10, 2),
            k=np.full(10, 0),
            thickness=0.05,
            coherence=False
        ),
        make_layer(
            n=np.full(10, 1),
            k=np.full(10, 0),
            thickness=0,
            coherence=False
        )
    ]

    output = SMMO(layers, config)()
    print(output)


if __name__ == "__main__":
    main()
