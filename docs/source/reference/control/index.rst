.. currentmodule:: cana.boolean_network

Control
=======

These are methods and modules used to calculate control in Boolean Networks. They are divided in dynamics- and structure-based methods.

Note that these methods do not need to be called directly, as :class:`cana.boolean_network.BooleanNetwork` provides the appropriate methods.
A typical control workflow starts by instantiating a :class:`.BooleanNetwork`, computing its attractors, and then evaluating which driver nodes are necessary to steer the system between those attractors. The reference below documents the individual helpers that make up that workflow.

.. code-block:: python

    from cana.boolean_network import BooleanNetwork

    # Define a two-node toggle switch where each gene inhibits the other.
    logic = {
        0: {"name": "A", "in": [1], "out": ["1", "0"]},
        1: {"name": "B", "in": [0], "out": ["1", "0"]},
    }

    bn = BooleanNetwork.from_dict(logic, name="Toggle switch")
    driver_sets = bn.attractor_driver_nodes()
    print(driver_sets)

The resulting ``driver_sets`` list contains the minimal node index combinations capable of reaching every attractor—each entry represents one feasible set of driver variables. See also the :doc:`../datasets/networks` catalogue for curated biological networks that you can explore with the same approach.

.. contents:: Contents
	:depth: 3

Dynamics based control
------------------------

The control methods used here are implemented directly on the base class :class:`.BooleanNetwork` and :class:`.BooleanNode`. That is because the Network class can ask its nodes directly to step into a specific trajectory, thus compartmentalizing the logic.

Attractor Control	
^^^^^^^^^^^^^^^^^

.. automethod:: cana.boolean_network.BooleanNetwork.attractor_driver_nodes
	:noindex:

.. automethod:: cana.boolean_network.BooleanNetwork.controlled_state_transition_graph
	:noindex:

.. automethod:: cana.boolean_network.BooleanNetwork.controlled_attractor_graph
	:noindex:

Structure based control
-----------------------

These are control methods that only take the structure of the boolean network (aka: the structure graph) into consideration when computing driver nodes.

.. automodule:: cana.control.fvs
	:members:
	
.. automodule:: cana.control.mds
	:members:

.. automodule:: cana.control.sc
	:members: