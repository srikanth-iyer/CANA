Drawing
=======

The :mod:`cana.drawing` module gathers utilities for visualising canalization
structures and Boolean node behaviour. These helpers rely on external
visualisation backends, so make sure the required dependencies are available
before calling them.

.. note::

   The canalizing map renderer depends on both the Graphviz **system
   package** (providing the ``dot`` command line tools) and the
   :mod:`graphviz` **Python package** that ships the bindings used by
   ``cana``. Install the operating system package through your package
   manager (for example ``apt install graphviz`` or ``brew install
   graphviz``) and then install the Python bindings with ``pip install
   graphviz``.

Canalizing maps
---------------

.. currentmodule:: cana.drawing.canalizing_map

.. autofunction:: draw_canalizing_map_graphviz

Look-up tables
--------------

.. currentmodule:: cana.drawing.plot_look_up_table

.. autofunction:: plot_look_up_table

.. autofunction:: plot_annigen_look_up_table
