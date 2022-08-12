"""
=====================
dna_sim_controller.py
=====================

This file defines the main controller for the simulation.
"""

from ..interface import Controller, DataContainer
from ..environment.dna_sim_environment import DNASimEnvironment
from ..variables import stage_per_collection, dt, total_time, scaling
from ..datacontainer.data_recorder import RNAPPositionRecorder


class RunConfig(DataContainer):
    def __init__(self, controller=None, record_rnap_position: bool = True, record_rnap_amount: bool = True,
                 record_rnap_state: bool = True, record_processing_time: bool = True,
                 record_protein_amount: bool = True, record_protein_production: bool = True,
                 log_data: bool = True, show_progress_bar: bool = True):
        super().__init__(parent=controller)
        self.record_rnap_position = record_rnap_position
        self.record_rnap_amount = record_rnap_amount
        self.record_rnap_state = record_rnap_state
        self.record_processing_time = record_processing_time
        self.record_protein_amount = record_protein_amount
        self.record_protein_production = record_protein_production
        self.log_data = log_data
        self.show_progress_bar = show_progress_bar


class DNASimController(Controller):
    """
    This is the central controller for the DNA simulation.
    """
    def __init__(self,  rnap_loading_rate: float, config: RunConfig = RunConfig(), **kwargs):
        super().__init__()
        # Setup
        self.time_index = 0
        self.env = DNASimEnvironment(controller=self, rnap_loading_rate=rnap_loading_rate, **kwargs)
        self.total_time = scaling(total_time)
        self.dt = dt
        self.stage_per_collection = stage_per_collection

        # Data Recording Setup
        self.config = config

        # Initialize Data Recorders
        self.position = RNAPPositionRecorder(self, self.env.dna, self.stage_per_collection, self.dt, self.total_time+1)
        pass

    def start(self):
        print("Init!")
        self.init()
        print("Start!")
        for i in range(self.total_time + 1):
            self.env.step(time_index=self.time_index)
            self._log()
            self.time_index += 1

        pass

    def init(self):
        self.env.init()
        pass

    def call_back(self, option, data):
        pass

    def _log(self):
        self.position.log(self.time_index)
        pass
