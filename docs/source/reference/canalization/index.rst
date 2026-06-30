.. currentmodule:: cana.canalization

Canalization
==============

Canalization tools in CANA reveal how Boolean network logic can be simplified by
identifying redundant or dispensable inputs. Use this module when you want to
quantify how stable a node's behavior is under perturbations, compare the
controllability of different subnetworks, or derive reduced representations that
highlight decisive control structures before running large-scale simulations.

The tutorials in the project repository walk through full analysis workflows—see
`Canalization - BioModels.ipynb <https://github.com/rionbr/CANA/blob/master/tutorials/Canalization%20-%20BioModels.ipynb>`_.
For a catalogue of functions and parameters refer to the
:mod:`cana.canalization.boolean_canalization` API reference.

.. rubric:: Quick start

.. code-block:: python

    from cana.datasets import bio

    # Load a bundled biological network and inspect the first node
    network = bio.MARQUESPITA()
    node = network.nodes[0]

    # Calculate a canalization measure (input redundancy) and print it
    redundancy = node.input_redundancy()
    print(f"{node.name} input redundancy: {redundancy:.3f}")

.. automodule:: cana.canalization.boolean_canalization
        :members:
