CANAlization: Control & Redundancy in Boolean Networks
=======================================================

This package implements a series of methods used to study control, canalization, and redundancy in Boolean Networks.

[![PyPI version](https://badge.fury.io/py/cana.svg)](https://badge.fury.io/py/cana)
[![Downloads](https://pepy.tech/badge/cana)](https://pepy.tech/project/cana)

## Citation

If you use `cana` in your research, please cite us and check out our related papers below!

- A.M. Marcus, J.C. Rozum, H. Sizek, L.M. Rocha [2025]. "[CANA v1.0.0: efficient quantification of canalization in automata networks](https://doi.org/10.1093/bioinformatics/btaf461)". *Bioinformatics*. btaf461. doi: 10.1093/bioinformatics/btaf461

## Installation

**Stable Release (from PyPI)**
```bash
pip install cana
```

**Development Release (from GitHub)**
```bash
pip install git+https://github.com/CASCI-lab/CANA
```

## Quick Start

Here's a quick example of how to create a `BooleanNetwork`, find its attractors, and compute its Dynamics Canalization Map (DCM).

```python
from cana.boolean_network import BooleanNetwork

# Define a simple 3-node network
logic = {
    0: {'name': 'A', 'in': [1, 2], 'out': [0, 1, 1, 0]},
    1: {'name': 'B', 'in': [0], 'out': [1, 0]},
    2: {'name': 'C', 'in': [0, 1], 'out': [0, 0, 1, 1]}
}

# Create the BooleanNetwork object
bn = BooleanNetwork.from_dict(logic, name='My Simple Network')

# Find attractors using the State Transition Graph (STG)
attractors = bn.attractors(mode='stg')
print(f'Attractors: {attractors}')

# Compute the Dynamics Canalization Map (DCM)
dcm = bn.dynamics_canalization_map()
print(f'DCM: {dcm}')
```

## Documentation

The full documentation, including detailed tutorials and API references, can be found at: [casci-lab.github.io/CANA/](https://casci-lab.github.io/CANA/)

## Development

`CANA` uses `Cython` to accelerate some of its computations. If you are contributing to the project and need to modify the `.pyx` files, you will need to have `Cython` installed in your development environment.

**Setting up a development environment:**
```bash
# Clone the repository
git clone https://github.com/CASCI-lab/CANA.git
cd CANA

# Install in editable mode with development dependencies
pip install -e .
```

The `setup.py` file is configured to automatically detect if `Cython` is available and will compile the `.pyx` files accordingly. If you modify any `.pyx` files, they will be re-compiled the next time you run the installation command.

**Running Tests:**
To run the test suite, use `pytest`:
```bash
pytest
```

## Credits

`CANA` was originally written by Rion Brattig Correia and Alexander Gates, and has been developed with the help of many others. Thanks to everyone who has improved `CANA` by contributing code, bug reports (and fixes), documentation, and input on design and features.

## Support

`CANA` has been developed with support from the [Complex Adaptive Systems and Computational Intelligence (CASCI)](https://homes.luddy.indiana.edu/rocha/casci.php) lab at Indiana University, Bloomington, and the [CAPES Foundation](https://www.gov.br/capes/pt-br) of the Ministry of Education of Brazil.