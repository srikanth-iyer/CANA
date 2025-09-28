import unittest
from cana.boolean_network import BooleanNetwork
from cana.datasets.bio import THALIANA

def test_EG_weight_THALIANA():
    """Test that effective graph in-degree edge weights are computed correctly."""
    network = THALIANA()
    network.effective_graph()

    true = []
    for i, node in enumerate(network.nodes):
        # get sum from nx object
        edgews = {edge: network._eg.edges[edge]["weight"] for edge in network._eg.edges if edge[1]==i}
        true.append(sum(edgews.values()))
    assert network.effective_indegrees() == sorted(true, reverse=True)

class TestBooleanNetwork(unittest.TestCase):
    def setUp(self):
        """Set up a simple Boolean Network for testing."""
        self.logic = {
            0: {'name': 'A', 'in': [1, 2], 'out': [0, 1, 1, 0]},
            1: {'name': 'B', 'in': [0], 'out': [1, 0]},
            2: {'name': 'C', 'in': [0, 1], 'out': [0, 0, 1, 1]}
        }
        self.bn = BooleanNetwork.from_dict(self.logic, name='Test Network')

    def test_initialization(self):
        """Test that the BooleanNetwork is initialized correctly."""
        self.assertEqual(self.bn.Nnodes, 3)
        self.assertEqual(self.bn.name, 'Test Network')
        self.assertEqual([node.name for node in self.bn.nodes], ['A', 'B', 'C'])

    def test_step(self):
        """Test a single step of the network."""
        self.assertEqual(self.bn.step('000'), '010')
        self.assertEqual(self.bn.step('110'), '101')

    def test_attractors_stg(self):
        """Test finding attractors using the STG method."""
        attractors = self.bn.attractors(mode='stg')
        self.assertEqual(len(attractors), 1)
        self.assertEqual(len(attractors[0]), 1)
        self.assertCountEqual(attractors[0], [5])

    def test_attractors_bns(self):
        """Test finding attractors using the BNS method."""
        attractors = self.bn.attractors(mode='bns')
        if attractors:
            self.assertEqual(len(attractors), 1)
            self.assertEqual(len(attractors[0]), 1)
            self.assertCountEqual(attractors[0], [5])

    def test_structural_graph(self):
        """Test the creation of the structural graph."""
        sg = self.bn.structural_graph()
        edges = list(sg.edges())
        self.assertCountEqual(edges, [(1, 0), (2, 0), (0, 1), (0, 2), (1, 2)])
