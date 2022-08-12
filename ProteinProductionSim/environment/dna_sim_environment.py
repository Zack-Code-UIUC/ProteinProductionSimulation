"""
===================
dna_sim_environment
===================
this file contains the Environment subclass that is used for dna_sim_controller.
"""

from ..interface import Environment
from ..entity.dna_strand import DNAStrand
from ..entity.rnap import RNAP


class DNASimEnvironment(Environment):
    """
    This class represents the environment which the DNA-to-protein Experiment is constructed.

    Parameters
    ----------


    Attributes
    ----------


    """
    def __init__(self, controller, rnap_loading_rate: float, **kwargs):
        super().__init__(parent=controller)
        self.dna = DNAStrand(environment=self, rnap_loading_rate=rnap_loading_rate, **kwargs)

        pass

    def step(self, time_index):
        return self.dna.step(time_index)

    def init(self):
        self.dna.init()
        pass

    def call_back(self, option, data):
        pass
