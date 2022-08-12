"""
===================
dna_sim_environment
===================
this file contains the Environment subclass that is used for dna_sim_controller.
"""

from ..interface import Environment
from ..entity.dna_strand import DNAStrand

class DNASimController(Environment):
    """
    This class represents the environment which the DNA-to-protein Experiment is constructed.
    """