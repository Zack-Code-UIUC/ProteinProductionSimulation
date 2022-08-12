"""
============
dna_strand.py
============

This file contains the DNAStrand class.
"""

from interface import Entity
from variables import length, scaling, dt, total_time, v_0, pauseDuration, pauseSite, RNAP_size, t_on, tau_0, tau_c

from helper.loading_list import LoadingList
from helper.supercoilling import s, n_dependence_cubic_3
from .rnap import RNAP
import numpy as np


class DNAStrand(Entity):
    """
    This class represents the DNA strand.

    Parameters
    ----------


    Attributes
    ----------
    """
    def __init__(self, environment, rnap_loading_rate, include_supercoiling=True, include_busty_promoter=False,
                 rnap_loading_pattern="stochastic", promoter_shut_off_time=-1, pause_profile="flat",
                 ribo_loading_profile="stochastic", degradation_profile="exponential",
                 protein_production_off: bool = False, if_storing_supercoiling_value: bool = False):
        super().__init__(environment)
        self.length = length
        self.rnap_loading_rate = rnap_loading_rate

        # settings
        self.include_supercoiling = include_supercoiling
        self.include_busty_promoter = include_busty_promoter
        self.rnap_loading_pattern = rnap_loading_pattern
        self.pause_profile = pause_profile
        self.ribo_loading_pattern = ribo_loading_profile
        self.degradation_profile = degradation_profile
        self.protein_production_off = protein_production_off

        # variables related to the status of its children mRNAs.
        self.loaded = 0  # number of RNAPs that have been loaded.
        self.detached = 0  # number of RNAPs that have been detached.
        self.attached = 0  # number of RNAPs that are attached.
        self.degrading = 0  # number of RNAPs that are degrading.
        self.degraded = 0  # number of RNAPs that are already degraded.

        # setup loading list.
        if_stochastic = False
        match self.rnap_loading_pattern:
            case "stochastic":
                if_stochastic = True
            case "uniform":
                if_stochastic = False
        if self.include_busty_promoter:
            self.loading_list = LoadingList(self, scaling(total_time), self.rnap_loading_rate*dt,
                                            if_stochastic=if_stochastic, if_bursty=True)
        else:
            self.loading_list = LoadingList(self, scaling(total_time), self.rnap_loading_rate*dt,
                                            if_stochastic=if_stochastic, if_bursty=False)

        # promoter_state
        self.promoter_state = False
        if promoter_shut_off_time == -1:
            self.T_stop = scaling(total_time)
        elif promoter_shut_off_time >= 0:
            self.T_stop = scaling(promoter_shut_off_time)

        # Site-Specific Pausing
        self.include_site_specific_pausing = True
        match self.pause_profile:
            case "flat":
                self.include_site_specific_pausing = False
            case "OnepauseAbs":
                self.n_pause = 1
            case "TwopauseAbs":
                self.n_pause = 2

        # mRNA Degradation
        self.RNAP_LIST = []

        # adaptive supercoiling
        self.loading_list.trim(self.T_stop)
        self.r_ref = np.zeros(self.loading_list.get_length())
        self.flag_r_ref = np.zeros(self.loading_list.get_length(), dtype=bool)
        self.T_open = scaling(total_time)
        self.just_loaded = False

        # storage
        self.if_storing_supercoiling_value = if_storing_supercoiling_value
        self.phi = None
        self.torq = None
        self.f_n = None
        self.velo = None
        self.stepping = None

    def init(self):
        for rnap in self.RNAP_LIST:
            rnap.init()
        pass

    def step(self, time_index):

        # REASON: time to check for loading
        to_load = False

        # REASON: check if there is one loading attempt, obviously, the after T_stop, nothing should load.
        if self.loading_list.if_can_load(time_index) and time_index <= self.T_stop:
            to_load = True

        # REASON: check if there is one RNAP congesting the loading site, we just check the last rnap.
        if to_load and len(self.RNAP_LIST) != 0:
            if (self.RNAP_LIST[-1].position - RNAP_size) < 0:
                to_load = False

        # REASON: if it can load, then load one RNAP
        if to_load:
            self.RNAP_LIST.append(
                RNAP(self, time_index, pause_profile=self.pause_profile, ribo_loading_profile=self.ribo_loading_pattern,
                     degradation_profile=self.degradation_profile, protein_production_off=self.protein_production_off))
            self.loaded += 1
            self.attached += 1
            if self.include_supercoiling:
                self.T_open = time_index + scaling(t_on)
                self.promoter_state = True
                self.just_loaded = True
                if self.loaded > 1:
                    self.r_ref[self.loaded - 2] = self.RNAP_LIST[self.loaded - 2].position
                    self.flag_r_ref[self.loaded - 2] = True

        # REASON: check the promoter closing resulting from the RNAP loading
        if self.just_loaded and time_index >= self.T_open and self.promoter_state:
            # REASON: switch the promoter
            self.just_loaded = False
            self.promoter_state = False
            if self.loaded != 0:
                # REASON: we reset the reference position of the last RNAP
                self.r_ref[self.loaded - 1] = self.RNAP_LIST[self.loaded - 1].position

        # REASON: check for permanent promoter shutoff
        if time_index >= self.T_stop and self.promoter_state:
            self.promoter_state = False
            if self.loaded != 0:
                value = self.RNAP_LIST[self.loaded - 1].position
                self.r_ref[self.loaded - 1] = value

        # REASON: calculate stepping
        if self.include_supercoiling:
            stepping = self.supercoiling()
        else:
            stepping = []
            for i in range(len(self.RNAP_LIST)):
                stepping.append(v_0*dt)

        # REASON: check for site-specific pausing and set pausing.
        if self.include_site_specific_pausing and len(self.RNAP_LIST) != 0:
            for count, rnap in enumerate(self.RNAP_LIST):
                if not rnap.attached:
                    # REASON: this just skip the detached RNAP. No need to check them
                    continue

                # REASON: if passed both, then just skip
                if rnap.passed_site_1 and rnap.passed_site_2:
                    continue

                if not rnap.passed_site_1:
                    # REASON: check if the rnap is still in the pause site.
                    if rnap.passing_1 and rnap.position + 1 / scaling(pauseDuration[0]) < pauseSite[0]:
                        stepping[count] = 1 / scaling(pauseDuration[0])
                        continue

                    # REASON: passing out of the pausing site.
                    if rnap.passing_1 and rnap.position + 1 / scaling(pauseDuration[0]) >= pauseSite[0]:
                        self.RNAP_LIST[count].passing_1 = False
                        self.RNAP_LIST[count].passed_site_1 = True
                        continue

                    # REASON: passing into the pause_site 1
                    if rnap.position < (pauseSite[0] - 1) and ((rnap.position + stepping[count]) >= pauseSite[0] - 1):
                        self.RNAP_LIST[count].passing_1 = True
                        self.RNAP_LIST[count].position = pauseSite[0] - 1
                        stepping[count] = 0
                        # print("Enter Pause Site 1")
                        continue

                if not rnap.passed_site_2:
                    if rnap.passing_2 and rnap.position + 1 / scaling(pauseDuration[1]) < pauseSite[1]:
                        stepping[count] = 1 / scaling(pauseDuration[1])
                        continue

                    # REASON: passing out of the pausing site.
                    if rnap.passing_2 and rnap.position + 1 / scaling(pauseDuration[1]) >= pauseSite[1]:
                        self.RNAP_LIST[count].passing_2 = False
                        self.RNAP_LIST[count].passed_site_2 = True
                        continue

                    # REASON: passing into the pause_site 2
                    if rnap.position < (pauseSite[1] - 1) and ((rnap.position + stepping[count]) >= pauseSite[1] - 1):
                        self.RNAP_LIST[count].passing_2 = True
                        self.RNAP_LIST[count].position = pauseSite[1] - 1
                        stepping[count] = 0
                        continue

        # REASON: check for hindrance and modify stepping
        for count, rnap in enumerate(self.RNAP_LIST):
            # REASON: if this rnap is not attached, we do not need to worry about it.
            if not rnap.attached:
                continue
            # REASON: if this rnap is the front-most attached rnap, we also do not need to care about it.
            if count == self.detached:
                continue
            # REASON: if this rnap is going to move past its front rnap, then a collision happens,
            #         and we need to set its position such that it sit right between its front rnap with a distance of
            #         RNAP_size.
            else:
                if (rnap.position + stepping[count]) > \
                        (self.RNAP_LIST[count - 1].position + stepping[count - 1] - RNAP_size):
                    stepping[count] = self.RNAP_LIST[count - 1].position + stepping[
                        count - 1] - RNAP_size - rnap.position
        # print(stepping)
        # REASON: now we plug in the stepping into the RNAPs and collect protein production from all RNAP.
        prot = 0

        for i in range(len(self.RNAP_LIST)):
            prot += self.RNAP_LIST[i].step(time_index, stepping[i])

        # this represents the amount of RNAPs that are attached at this moment.
        # print([self.loaded, self.attached, self.detached, self.degrading, self.degraded])

        # return protein production
        return prot

    def call_back(self, option: str, data=1):
        match option:
            case "attached":
                self.attached += data
            case "detached":
                self.detached += data
            case "degrading":
                self.degrading += data
            case "degraded":
                self.degraded += data

    def supercoiling(self):

        # STEP: setup
        size = self.attached

        # STEP: phi generation
        PHI = np.zeros(size)
        '''
        for i in range(size + 1):
            j = i + self.detached  # STEP: convert to real index
            if i == 0:  # STEP: check the front-most element
                PHI[i] = 0.0
            elif i == size:  # STEP: check the back-most element
                if self.promoter_state:  # STEP: check if the repressor is off. The supercoiling is diffused
                    PHI[i] = 0.0
                else:  # STEP: check if the repressor is attached, full supercoiling.
                    travel_dist_last_rnap = self.RNAP_LIST[j - 1].position - self.r_ref[j - 1]
                    PHI[i] = s(travel_dist_last_rnap)
            else:  # STEP: check the middle element
                travel_dist_last_rnap = self.RNAP_LIST[j - 1].position - self.r_ref[j - 1]
                PHI[i] = s(travel_dist_last_rnap - self.RNAP_LIST[j].position)
                # not passed the location
        '''
        for i in range(size):
            j = i + self.detached  # STEP: convert to real index
            positive = 0
            negative = 0
            R = self.RNAP_LIST[j].position - self.r_ref[j]
            if j == self.detached:  # STEP: check the front-most RNAP
                positive = 0
                if j == self.loaded - 1:
                    negative = s(R)
                else:
                    negative = s(R-self.RNAP_LIST[j+1].position)
            elif j == self.loaded - 1:  # STEP: check the back-most RNAP
                R_last = self.RNAP_LIST[j - 1].position - self.r_ref[j - 1]
                positive = s(self.RNAP_LIST[j].position - R_last)
                negative = s(R)
            else:  # STEP: check the middle element
                R_last = self.RNAP_LIST[j - 1].position - self.r_ref[j - 1]
                positive = s(self.RNAP_LIST[j].position - R_last)
                # print([i, j, self.loaded, self.detached])
                negative = s(R - self.RNAP_LIST[j+1].position)

            if self.promoter_state:
                negative = 0
            PHI[i] = negative + positive

        # STEP: torque generation
        torq = np.zeros(size)
        n = n_dependence_cubic_3(size)

        for i in range(size):
            torq[i] = tau_0 * n * (PHI[i])
        # STEP: velocity generation
        velo = np.zeros(len(self.RNAP_LIST))

        for i in range(self.detached):
            velo[i] = 0

        for i in range(
                size):  # we are using this for loop, to avoid large value generated from the exponential function.
            j = i + self.detached
            if torq[i] > 1.5 * tau_c:
                velo[j] = 0
            elif torq[i] < -1.5 * tau_c:
                velo[j] = 2 * v_0
            else:
                velo[j] = 2 * v_0 / (1 + np.exp(2 * (torq[i] / tau_c) ** 3))
        stepping = velo*dt
        if self.if_storing_supercoiling_value:
            self.phi = PHI
            self.torq = torq
            self.f_n = n
            self.velo = velo
            self.stepping = stepping
        return stepping
