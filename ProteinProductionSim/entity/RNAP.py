# RNAP Class

class RNAP:
    def __init__(self, t, pause_profile="flat", mRNA_degradation=True, DNA=None):
        # print("loaded")
        self.initial_t = t  # means the initial time when the RNAP attached
        self.position = 0  # means position =)
        self.initiated = False  # initiated indicates if the length has passed the size required for initiation (33nts)
        self.attached = True  # indicated if the RNAP is attached to the DNA. Will detached if reach the end
        self.dna = DNA  # this store the reference to its mother DNA, so that callback method can be used.

        # site-pausing
        self.passed_site_1 = False  # passed_site_1 is indicateing whether the RNAP has passed the pausing site. Also used to simply by-pass
        self.passed_site_2 = False
        self.passing_1 = False  # if the RNAP is passing 1 or 2.
        self.passing_2 = False
        if pause_profile == pauseProfile[0]:
            self.passed_site_1 = True
            self.passed_site_2 = True
        elif pause_profile == pauseProfile[1]:
            self.passed_site_1 = rand.choice(a=[True, False], p=[1 - pauseProb, pauseProb])
            self.passed_site_2 = True
        elif pause_profile == pauseProfile[2]:
            self.passed_site_1 = rand.choice(a=[True, False], p=[1 - pauseProb, pauseProb])
            self.passed_site_2 = rand.choice(a=[True, False], p=[1 - pauseProb, pauseProb])

        # mRNA degradation
        self.mRNA_degradation = mRNA_degradation
        self.degraded_length = 0
        self.degrading = False
        self.degraded = False
        if mRNA_degradation:
            self.t_degrade = rand.exponential(scale=ribo_loading_interval)
        else:
            self.t_degrade = total_time
        self.last_ribo_time_degraded = total_time

        # Loading of Ribosomes
        self.loading_list = RIBO_loading_list(self.t_degrade, kRiboLoading)
        if protein_production_off:
            self.loading_list = []
        self.loading_stage = 0  # tracks which index of loading list to check next.
        self.loading_number = len(self.loading_list)  # how many loading events are there in the loading_list.
        self.RIBO_LIST = np.zeros(self.loading_number)  # contains positions of all Ribosomes.
        self.RIBO_loaded = 0  # shows how many ribosomes are loaded onto the mRNA
        self.RIBO_detached = 0  # shows how many ribomes has detached.
        self.detach_time = total_time

    # RNAP.step() take in the base pair this instance transverse and change its inner data accordingly.
    def step(self, t, pace):

        if self.degraded:
            return 0
        prot = 0

        # check if detached
        if self.attached and (self.position + pace >= length):
            self.attached = False
            self.position = length + 1000
            pace = 0
            self.detach_time = t
            self.dna.detached += 1

        # recalculate(update) degradation
        if not self.degrading:  # if the mRNA is not degrading
            if t >= self.t_degrade + self.initial_t:  # if the time has exceed the degradation time, then the degrading state will be mark True.
                self.degrading = True
                self.dna.degrading += 1

            if not self.initiated:  # if not marked initiated, then check for if the mRNA has exceed that state. If so, mark the initiated True.
                if self.position + pace >= initiation_nt:
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

        # return protein production
        return prot

    # return internal data
    def get_internal(self):
        pass

# ||---------------------------------------------------------------------||

