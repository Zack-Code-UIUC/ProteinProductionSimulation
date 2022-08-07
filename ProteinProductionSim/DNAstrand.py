# ||---------------------------------------------------------------------||
# DNAstrand Class

# DNAstrand class handle the supercoiling, RNAP site-specific pausing and RNAP-to-RNAP hindrance and assign degradation coefficient to mRNA instance.
# RNAP class handle the ribosome progression and act as a data container for RNAP-specific stuffs.
class DNAstrand:

    def __init__(self, interval_loading, pause_profile="flat", degradation_time=0, include_supercoiling=True,
                 include_busty_promoter=True, mRNA_degradation=True):
        self.length = length
        self.total_size_loading = 0  # number of RNAPs loading events what will be checked. will be updated by the initiation()

        self.loaded = 0  # number of RNAPs that have been loaded.
        self.detached = 0  # number of RNAPs that have been detached.
        self.attached = 0  # number of RNAPs that are attached.
        self.degrading = 0  # number of RNAPs that are degrading.
        self.degraded = 0  # number of RNAPs that are already degraded.

        # Supercoiling
        self.include_supercoiling = include_supercoiling

        # Busty Promoter
        if include_busty_promoter:
            tau_on = tau_off * tau_loading / (interval_loading - tau_loading)
            self.loading_list = RNAP_loading_list(promoter_time_list(total_time, tau_on, tau_off), 1 / tau_loading)
        else:
            self.loading_list = RNAP_loading_list([[True, total_time]], 1 / interval_loading)
            self.promoter_list = []
        self.promoter_state = False
        self.T_stop = 90  # total_time

        # Site-Specific Pausing
        self.pause_profile = pause_profile
        self.include_site_specific_pausing = True
        if pause_profile == pauseProfile[0]:
            self.include_site_specific_pausing = False
            self.n_pause = 0
        elif pause_profile == pauseProfile[1]:
            self.n_pause = 1
        elif pause_profile == pauseProfile[2]:
            self.n_pause = 2

        # mRNA Degradation
        self.mRNA_degradation = mRNA_degradation
        self.t_degrade = degradation_time
        self.RNAP_LIST = []

        # Test
        self.print_once = True

    def initiation(self):
        # trim the rnap loading list such that after T stop, no loading happens.
        for i in self.RNAP_LIST:
            if i > self.T_stop:
                self.RNAP_LIST.remove(i)

        # recording
        if group_two:
            self.phi_first_rnap = np.zeros(int(total_time / dt))
            self.torq_first_rnap = np.zeros(int(total_time / dt))
            self.f_n = np.zeros(int(total_time / dt))
            self.record_promoter_state = np.zeros(int(total_time / dt))
            self.record_last_r_ref = np.zeros(int(total_time / dt))

        self.total_size_loading = len(self.loading_list)
        self.promoter_state = False

        # Adaptive Supercoiling
        self.r_ref = np.zeros(self.total_size_loading)
        self.flag_r_ref = np.zeros(self.total_size_loading, dtype=bool)
        self.T_open = total_time
        self.just_loaded = False

    def step(self, t, index):
        self.index = index

        # check the promoter closing resulting from the RNAP loading
        if self.just_loaded and t > self.T_open and self.promoter_state:
            self.just_loaded = False
            self.promoter_state = False
            if self.loaded != 0:
                value = self.RNAP_LIST[self.loaded - 1].position
                self.r_ref[self.loaded - 1] = value

        # check for promoter state
        if t + dt > self.T_stop and self.promoter_state:
            self.promoter_state = False
            if self.loaded != 0:
                value = self.RNAP_LIST[self.loaded - 1].position
                self.r_ref[self.loaded - 1] = value

        # calculate supercoiling
        if self.include_supercoiling:
            if adaptive_supercoiling:
                velo = self.adaptive_supercoiling()
            else:
                velo = self.supercoiling()
        else:
            velo = []
            for i in range(len(self.RNAP_LIST)):
                velo.append(v_0)

        # check for site-specific pausing and set pausing.
        if self.include_site_specific_pausing and len(self.RNAP_LIST) != 0:
            for count, rnap in enumerate(self.RNAP_LIST):
                if rnap.passed_site_1 and rnap.passed_site_2:
                    continue

                if not rnap.passed_site_1:
                    # check if the the rnap is still in the pause site.
                    if rnap.passing_1 and rnap.position + 1 / pauseDuration[0] < pauseSite[0]:
                        velo[count] = 1 / pauseDuration[0]
                        continue

                    # passing out of the pausing site.
                    if rnap.passing_1 and rnap.position + 1 / pauseDuration[0] >= pauseSite[0]:
                        self.RNAP_LIST[count].passing_1 = False
                        self.RNAP_LIST[count].passed_site_1 = True
                        continue

                    # passing into the pause_site 1
                    if rnap.position < (pauseSite[0] - 1) and ((rnap.position + dt * velo[count]) >= pauseSite[0] - 1):
                        self.RNAP_LIST[count].passing_1 = True
                        self.RNAP_LIST[count].position = pauseSite[0] - 1
                        velo[count] = 0
                        # print("Enter Pause Site 1")
                        continue

                if not rnap.passed_site_2:
                    if rnap.passing_2 and rnap.position + 1 / pauseDuration[1] < pauseSite[1]:
                        velo[count] = 1 / pauseDuration[1]
                        continue

                    # passing out of the pausing site.
                    if rnap.passing_2 and rnap.position + 1 / pauseDuration[1] >= pauseSite[1]:
                        self.RNAP_LIST[count].passing_2 = False
                        self.RNAP_LIST[count].passed_site_2 = True
                        continue

                    # passing into the pause_site 2
                    if rnap.position < (pauseSite[1] - 1) and ((rnap.position + dt * velo[count]) >= pauseSite[1] - 1):
                        self.RNAP_LIST[count].passing_2 = True
                        self.RNAP_LIST[count].position = pauseSite[1] - 1
                        velo[count] = 0
                        continue

        # calculate the stepping
        stepping = []
        for velocity in velo:
            stepping.append(velocity * dt)

        # check for hindrance and modify stepping
        for count, rnap in enumerate(self.RNAP_LIST):
            if count == 0:
                continue
            if not rnap.attached:
                continue
            else:
                if (rnap.position + stepping[count]) > (
                        self.RNAP_LIST[count - 1].position + stepping[count - 1] - RNAP_size):
                    stepping[count] = self.RNAP_LIST[count - 1].position + stepping[
                        count - 1] - RNAP_size - rnap.position

        # plug the steppings into the RNAPs and collect protein production from all RNAP.
        # implement Multiprocessing here!
        # 6 processes at once in max.

        prot = 0

        for i in range(len(self.RNAP_LIST)):
            prot += self.RNAP_LIST[i].step(t, stepping[i])
        self.attached = self.loaded - self.detached  # this represent the amount of rnaps that are attached at this moment.

        # time to check for loading
        to_load = False

        # check if there is one loading attempt
        if len(self.loading_list) != 0 and t + dt >= self.loading_list[0] and t + dt <= self.T_stop:
            to_load = True

        # check if there is one RNAP congesting the loading site
        if to_load and len(self.RNAP_LIST) != 0:
            if (self.RNAP_LIST[-1].position - RNAP_size) < 0:
                to_load = False

        # if can load, then load one RNAP
        if to_load:
            self.RNAP_LIST.append(
                RNAP(t + dt, pause_profile=self.pause_profile, mRNA_degradation=self.mRNA_degradation, DNA=self))
            self.loading_list.pop(0)
            self.loaded = len(self.RNAP_LIST)
            if adaptive_supercoiling:
                self.T_open = t + t_on
                self.promoter_state = True
                self.just_loaded = True
                if self.loaded > 1:
                    self.r_ref[self.loaded - 2] = self.RNAP_LIST[self.loaded - 2].position
                    self.flag_r_ref[self.loaded - 2] = True
        self.attached = self.loaded - self.detached  # this represent the amount of rnaps that are attached at this moment.

        # record
        if group_two:
            self.record_promoter_state[index] = self.promoter_state
            if self.loaded != 0:
                self.record_last_r_ref[0] = self.r_ref[0]

        # return protein production
        return prot

    def supercoiling(self):
        # setup
        size = self.attached

        # phi generation
        PHI = np.zeros(size + 1)

        for i in range(size + 1):
            j = i + self.detached
            if i == 0:
                PHI[i] = S(self.RNAP_LIST[j].position)
            elif i == size:
                PHI[i] = S(self.RNAP_LIST[j - 1].position)
            else:
                PHI[i] = S(self.RNAP_LIST[j - 1].position - self.RNAP_LIST[j].position)

        # torque generation
        torq = np.zeros(size)
        n = n_dependence_cubic_3(size)
        for i in range(size):
            torq[i] = -tau_0 * n * (PHI[i] - PHI[i + 1])

        # velocity generation
        size = len(self.RNAP_LIST)
        velo = np.zeros(size)

        for i in range(self.detached):
            velo[i] = 0

        for i in range(np.size(torq)):
            j = i + self.detached
            if torq[i] > 1.5 * tau_c:
                velo[j] = 0
            elif torq[i] < -1.5 * tau_c:
                velo[j] = 2 * v_0
            else:
                velo[j] = 2 * v_0 / (1 + np.exp(2 * (torq[i] / tau_c) ** 3))
        if group_two:
            self.phi_first_rnap[self.index] = PHI[1]
            if self.RNAP_LIST[0].attached:
                self.torq_first_rnap[self.index] = torq[0]
            self.f_n[self.index] = n
        return velo

    def adaptive_supercoiling(self):
        ''''''
        # setup
        size = self.attached

        # phi generation
        PHI = np.zeros(size + 1)

        for i in range(size + 1):
            j = i + self.detached  # real index
            if i == 0:  # frontmost element
                PHI[i] = 0.0
            elif i == size:  # backmost element
                if self.promoter_state:  # the repressor is off. The supercoiling is diffused
                    PHI[i] = 0.0
                else:  # the repressor is attached, full supercoiling.
                    travel_dist_last_rnap = self.RNAP_LIST[j - 1].position - self.r_ref[j - 1]
                    PHI[i] = S(travel_dist_last_rnap)
            else:  # middle element
                travel_dist_last_rnap = self.RNAP_LIST[j - 1].position - self.r_ref[j - 1]
                PHI[i] = S(travel_dist_last_rnap - self.RNAP_LIST[j].position)
                # not passed the location

        # torque generation
        torq = np.zeros(size)
        n = n_dependence_cubic_3(size)
        for i in range(size):
            torq[i] = -tau_0 * n * (PHI[i] - PHI[i + 1])

        # velocity generation
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

        # record
        if group_two:
            if self.loaded != 0 and self.RNAP_LIST[0].attached:
                self.torq_first_rnap[self.index] = torq[0]
                self.phi_first_rnap[self.index] = PHI[1]
            self.f_n[self.index] = n
        return velo

# ||---------------------------------------------------------------------||
