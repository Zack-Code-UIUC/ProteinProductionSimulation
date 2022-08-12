"""
================
data_recorder.py
================
"""
import numpy as np

from interface import DataContainer, Environment, Controller
import matplotlib.pyplot as plt
from entity.dna_strand import DNAStrand
from environment.dna_sim_environment import DNASimEnvironment


class DataRecorder(DataContainer):
    """
    This is the basis class for the data recorder.
    """
    def __init__(self, controller: Controller, collection_interval):
        super().__init__(controller)
        self.collection_interval = collection_interval
        pass

    def init(self):
        pass

    def log(self, time_index):
        pass

    def get(self):
        pass

    def plot(self, fig):
        pass

    def store(self, path):
        pass


class RNAPPositionRecorder(DataRecorder):
    def __init__(self, controller, target: DNAStrand, collection_interval: int, dt: float, total_time: int):
        super().__init__(controller, collection_interval)
        self._tot_time = total_time
        self._dna = target
        self._dt = dt
        self._data = []
        self._length = 0
        self._rnap_loaded = 0
        self._rnap_detached = 0
        self._start_end = []

    def log(self, time_index: int):
        # STEP: check RNAP amount through loaded and attached
        if time_index % self.collection_interval == 0:
            # STEP: check if detached increase
            if self._dna.detached > self._rnap_detached:
                self._start_end[self._rnap_detached][1] = time_index - 1
                self._rnap_detached += 1
            # STEP: check if loaded increase
            if self._dna.loaded > self._rnap_loaded:
                self._start_end.append([time_index, -1])
                self._rnap_loaded += 1
                self._data.append(np.zeros(self._tot_time))

        # STEP: now we track the position.
        for i in range(self._rnap_loaded):
            if i < self._rnap_detached:
                continue
            self._data[i][self._length] = self._dna.RNAP_LIST[i].position

        self._length += 1
        pass

    def plot(self, axe):
        axe.set_xlabel('Time [s]')
        axe.set_ylabel('Position [bps]')
        axe.set_title('RNAP Position Plot')
        axe.grid(True)
        for i in range(self._rnap_loaded):
            start_index = self._start_end[i][0]
            stop_index = self._start_end[i][1]
            if stop_index < 0:
                stop_index = self._tot_time
            time_list = np.arange(start=start_index, stop=stop_index,
                                  step=1, dtype=float)
            axe.plot(time_list*self._dt, self._data[i][start_index:stop_index])
        print([self._rnap_loaded, self._rnap_detached])
        pass


# Amount Recorder Single-Value Recorder.

# Multi-Value Recorder
