"""
=======
mrna.py
=======

This file represent the messanger RNA, abbreviated as mRNA. The mRNA is produced from the transcription of DNA by the
RNAP. Ribosomes will attach to mRNA to conduct translation to synthesize proteins.


"""
import numpy as np
from ..interface import Entity
from ..helper.random_generator import exponential_generator, stepwise_exponential_generator
from ..helper.loading_list import LoadingList
from ..variables import kRiboLoading, ribo_loading_interval, protein_production_off, m1, m2, t_crit, initiation_nt, \
    RIBO_size


class mRNA(Entity):
    """
    Represent one mRNA. The mRNA class also manage all the ribosome and protein related action.
    """
    loadingProfile = ["uniform", "stochastic"]
    degradationProfile = ["exponential", "stepwise exponential"]

    def __init__(self, rnap, loading_profile="stochastic", degradation_profile="exponential"):
        self.parent = rnap  # this is the parent RNAP.
        self.length = 0  # this is the length of the mRNA
        self.initiated = False  # initiated indicates if the length has passed the size required for initiation (33nts)

        # mRNA degradation
        self.mRNA_degradation = degradation_profile
        self.degrading = False
        self.degraded = False
        match self.mRNA_degradation:
            case "exponential":
                self.t_degrade = exponential_generator(1/ribo_loading_interval)
            case "stepwise exponential":
                self.t_degrade = stepwise_exponential_generator(m1, m2, t_crit)

        # Loading of Ribosomes
        match loading_profile:
            case "uniform":
                self.loading_list = LoadingList(self, self.t_degrade, kRiboLoading, if_stochastic=False, if_bursty=False)
            case "stochastic":
                self.loading_list = LoadingList(self, self.t_degrade, kRiboLoading, if_stochastic=True, if_bursty=False)

        if protein_production_off:
            self.loading_list.dump()
        self.loading_stage = 0  # tracks which index of loading list to check next.
        self.loading_number = len(self.loading_list)  # how many loading events are there in the loading_list.
        self.RIBO_LIST = np.zeros(self.loading_number)  # contains positions of all Ribosomes.
        self.RIBO_loaded = 0  # shows how many ribosomes are loaded onto the mRNA
        self.RIBO_detached = 0  # shows how many ribosomes has detached.
        pass

    def _init(self):
        pass

    def step(self, time):
        if self.degraded:
            return
        prot = 0

        # recalculate(update) degradation
        if not self.degrading:  # if the mRNA is not degrading
            if time >= self.t_degrade + self.initial_t:  # if the time has exceeded the degradation time,
                self.degrading = True                    # then the degrading state will be mark True.
                self.dna.degrading += 1

            if not self.initiated:  # if not marked initiated, then check for if the mRNA has exceeded that state.
                if self.position + pace >= initiation_nt:  # If so, mark the initiated True.
                    self.initiated = True

            # if the mRNA is initiated, then check for potential Loading.
            if self.loading_stage != self.loading_number and t >= self.initial_t + self.loading_list[
                self.loading_stage]:  # first check if the loading list has been exhausted, then check for potential loading.
                self.loading_stage += 1  # increment the loading state counter.
                if self.initiated:
                    if self.RIBO_loaded == 0:  # if there are zero Ribosome loaded, then no worry for hindrance, just load it.
                        self.RIBO_loaded += 1
                    else:
                        if self.RIBO_LIST[self.RIBO_loaded - 1] - RIBO_size > 0:  # check for hindrance.
                            self.RIBO_loaded += 1

        # set up stepping for Ribosome and check for hindrance
        active_RIBO = self.RIBO_loaded - self.RIBO_detached
        stepping = np.zeros(active_RIBO)
        for i in range(active_RIBO):
            stepping[i] = v_0 * dt

            # check for ribo hindrance and the potential ribosome crowding near the end if the mRNA is not detached.
        for i in range(active_RIBO):
            actual_i = i + self.RIBO_detached
            # check for ribo hindrance, with two potential case of hindrance:
            # 1. if the mRNA has not finished transcription, then the leading Ribosome cannot move ahead of that length transcribed, as that would be unphysical
            if i == 0 and self.attached:
                if self.RIBO_LIST[actual_i] + stepping[i] > self.position:
                    stepping[i] = self.position - self.RIBO_LIST[actual_i]
            # 2. hindrance between the adjacent Ribosomes
            if i != 0:
                if self.RIBO_LIST[actual_i] + stepping[i] > self.RIBO_LIST[actual_i - 1] + stepping[i - 1] - RIBO_size:
                    stepping[i] = self.RIBO_LIST[actual_i - 1] + stepping[i - 1] - RIBO_size - self.RIBO_LIST[actual_i]

        # update ribosome position and check for protein production and detached Ribosome
        old_RIBO_detached = self.RIBO_detached
        for i in range(active_RIBO):
            actual_i = i + old_RIBO_detached
            self.RIBO_LIST[actual_i] += stepping[i]
            if self.RIBO_LIST[actual_i] > length:
                self.RIBO_detached += 1
                prot += 1

        # now add pace
        self.position += pace

        # check for complete degradation.
        if not self.attached and self.degrading and not self.degraded and (self.RIBO_loaded == self.RIBO_detached):
            self.degraded = True
            self.dna.degraded += 1
            self.last_ribo_time_degraded = t

        pass

    def _call_back(self):
        pass


