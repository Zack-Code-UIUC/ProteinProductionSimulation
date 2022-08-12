"""
=================
ribo_container.py
=================

This file defines the RIBOContainer class which is used by the RNAP class to store information related to the ribosomes.

"""

from interface import DataContainer
from variables import RIBO_size, dt, k_elong, length
import numpy as np


class RIBOContainer(DataContainer):
    def __init__(self, rnap):
        super().__init__(rnap)
        self.ribo_loaded = 0
        self.ribo_detached = 0
        self.ribo_list = []

    def if_empty(self):
        if self.ribo_loaded == 0:
            return True
        return False

    def if_all_detached(self):
        return self.ribo_loaded == self.ribo_detached

    def if_clear_at_start(self):
        if self.if_empty():
            return True
        else:
            if self.ribo_list[self.ribo_loaded-1] - RIBO_size >= 0:
                return True
            else:
                return False

    def load_one(self):
        self.ribo_loaded += 1
        self.ribo_list.append(0)

    def step(self, time_index, rnap_position):

        # REASON: we now will create a general stepping array, that just represents the stepping for each attached
        #         ribosome. At first, the stepping for any attached ribosome is just the elongation speed times the dt,
        #         but we also have to account for hindrance between each ribosome.

        # REASON: set up the general stepping array.
        stepping = np.full(shape=self.ribo_loaded, fill_value=k_elong*dt, dtype=float)

        # REASON: check for hindrance. start from the front ribosome and adjust as we're moving toward the end.
        for i in range(self.ribo_loaded):
            # REASON: this is the actual index of the element within the self.ribo_list.
            actual_i = i + self.ribo_detached

            # REASON: check for ribo hindrance. there are two potential cases of hindrance:
            #         First: if the mRNA has not finished transcription, then the front-most Ribosome cannot move ahead
            #         of that length transcribed, as that would be unphysical
            if i == 0 and self.parent.attached:
                if self.ribo_list[actual_i] + stepping[i] > rnap_position:
                    stepping[i] = rnap_position - self.ribo_list[actual_i]
            # 2. hindrance between the adjacent Ribosomes
            if i != 0:
                if self.ribo_list[actual_i]+stepping[i] > self.ribo_list[actual_i - 1]+stepping[i - 1]-RIBO_size:
                    stepping[i] = self.ribo_list[actual_i-1]+stepping[i-1]-RIBO_size-self.ribo_list[actual_i]

        # update ribosome position and check for protein production and detached Ribosome
        prot = 0
        for i in range(self.ribo_loaded):
            actual_i = i + self.ribo_detached
            self.ribo_list[actual_i] += stepping[i]
            if self.ribo_list[actual_i] > length:
                self.ribo_detached += 1
                prot += 1

        return prot
