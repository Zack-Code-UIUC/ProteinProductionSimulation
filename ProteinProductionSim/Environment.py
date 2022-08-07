# ||---------------------------------------------------------------------||
# Environment Class

# The class Environment is used to encapsulate all activities, including initiation of the whole sytem.
# Responsible for stepping of all variable and class, including t.
# Also responsible for recording all the data and their output.
class Environment:

    def __init__(self, interval_loading=500, pause_profile="flat", degradation_time=20, include_supercoiling=True,
                 include_busty_promoter=True, mRNA_degradation=True):
        # record the profile
        self.interval_loading = interval_loading
        self.pause_profile = pause_profile
        self.include_supercoiling = include_supercoiling
        self.include_busty_promoter = include_busty_promoter
        self.mRNA_degradation = mRNA_degradation

        # setup the profile
        self.t = 0
        self.DNA = DNAstrand(self.interval_loading, pause_profile=self.pause_profile, degradation_time=20,
                             include_supercoiling=self.include_supercoiling,
                             include_busty_promoter=self.include_busty_promoter, mRNA_degradation=self.mRNA_degradation)
        self.total_prot = 0

    def step(self, index):  # Used by start().
        return self.DNA.step(self.t, index)

    def start(self):
        stage_per_collection = int(data_collection_interval / dt)
        self.initialize_data_structure()
        if testing: self.print_initialization()

        total_iteration = int(total_time / data_collection_interval)
        start_time = time.perf_counter()
        per_stage = int(data_collection_interval / dt)
        for i in range(int(total_time / data_collection_interval)):
            t_i = time.perf_counter()
            for j in range(stage_per_collection):
                prot = self.step(i * per_stage + j)
                # record protein production
                if group_one:
                    self.protein_production[i * per_stage + j] = prot
                self.total_prot += prot
                # progress time and stage
                self.t += dt
            self.step_record(i)
            t_f = time.perf_counter()
            # recording processing for each stage
            if group_four:
                self.step_processing_time[i] = t_f - t_i
            if testing: printProgressBar_RNAP_prot(i + 1, total_iteration, len(self.DNA.RNAP_LIST),
                                                   self.total_prot)  # print out an epicly cool real-time progress bar (^v^)
        end_time = time.perf_counter()
        self.total_runtime = end_time - start_time
        if testing: self.print_finish()
        pass

    # This helper function is created separatedly, because the user might want to do modification prior to using start(). By calling initialize_data_structure() inside the start(), the inconsistency from the modification is avoided.
    def initialize_data_structure(self):
        if group_one:
            # Protein production amount each stepping
            self.protein_cumulative_amount = np.zeros(int(total_time / data_collection_interval))

            # Protein production amount each stepping
            self.protein_production = np.zeros(int(total_time / dt))

            # RNAP amount
            self.RNAP_amount = np.zeros(int(total_time / data_collection_interval))

            # RNAP position
            self.RNAP_position = np.zeros((int(total_time / data_collection_interval), self.DNA.total_size_loading))

        if group_two:
            pass

        if group_three:
            # The amount of RNAPs in various states.
            self.rnap_loaded = np.zeros(int(total_time / data_collection_interval))
            self.rnap_detached = np.zeros(int(total_time / data_collection_interval))
            self.rnap_degrading = np.zeros(int(total_time / data_collection_interval))
            self.rnap_degraded = np.zeros(int(total_time / data_collection_interval))

        if group_four:
            # Processing time of each stage recorded
            self.step_processing_time = np.zeros(int(total_time / data_collection_interval))

        # DNA recording
        self.DNA_loading_list = np.copy(self.DNA.loading_list)

    def print_initialization(self):
        if not testing:
            return None
        initialization_bar_size = 80
        include_busty_promoter = True
        bar = '||' + '=' * int(
            (initialization_bar_size - 21) / 2 - 4) + '||' + 'Environment Parameters' + '||' + '=' * int(
            (initialization_bar_size - 21) / 2 - 4) + '||'
        print()
        print(bar)
        print('|  ', "{:<30} {:<15} {:<26}".format('name', 'value', ''), ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('interval_loading', str(self.interval_loading), ''), ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('include_supercoiling', str(self.include_supercoiling), ''), ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('include_busty_promoter', str(self.include_busty_promoter), ''),
              ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('include_mRNA_degradation', str(self.mRNA_degradation), ''), ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('length_rnap_loading_list', str(len(self.DNA_loading_list)), ''),
              ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('total simulation time', str(total_time) + 's', ''), ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('dt', "%.2f" % dt + 's', ''), ' |')
        print('|  ',
              "{:<30} {:<15} {:<26}".format('data_collection_interval', "%.2f" % data_collection_interval + 's', ''),
              ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('data_smoothing_interval', str(data_smoothing_interval) + 's', ''),
              ' |')
        print('|  ',
              "{:<30} {:<15} {:<26}".format('total collections', str(int(total_time / data_collection_interval)), ''),
              ' |')
        print('|  ', "{:<30} {:<15} {:<26}".format('stage per collection', str(int(data_collection_interval / dt)), ''),
              ' |')
        bar = '||' + '=' * (initialization_bar_size - 4) + '||'
        print(bar)
        print('Simulatation Starts. \n')

    def print_finish(self):
        if not testing:
            return None
        print(f'Simulation Runtime: {self.total_runtime} s.')
        # Simulation Preliminary Report

    # def dt_record(self, index):

    # Check and record data every step time as defined in the parameter. This step time is different from dt.
    def step_record(self, index):

        if group_one:

            # collect the protein amount
            self.protein_cumulative_amount[index] = self.total_prot
            # collect the rnap amount
            self.RNAP_amount[index] = self.DNA.attached

            # collect the rnap positions
            for i in range(len(self.DNA.RNAP_LIST)):
                self.RNAP_position[index][i] = self.DNA.RNAP_LIST[i].position

        if group_three:
            # collect the amount of RNAPs in various states.
            self.rnap_loaded[index] = self.DNA.loaded
            self.rnap_detached[index] = self.DNA.detached
            self.rnap_degrading[index] = self.DNA.degrading
            self.rnap_degraded[index] = self.DNA.degraded

    # Used for obtaining the inner store data. No need for this now, as all data are exposed.
    def get_internal(self):
        pass

    # Find the start index and end index of each RNAP
    # Data is outputed in the format of ((start_index, end_index, if_has_loaded, if_has_detached), (..., ..., ..., ...) ...)
    # Used as helper function by other functions.
    def find_start_end(self):
        start_end = np.zeros((self.RNAP_position.shape[1], 4))

        for i in range(self.RNAP_position.shape[1]):
            start = 0
            end = self.RNAP_position.shape[0]
            has_started = False
            has_ended = False
            for j in range(self.RNAP_position.shape[0]):
                if not has_ended:  # the RNAP has not yet detacheded
                    if not has_started:  # the RNAP has not attached
                        if not self.RNAP_position[j, i] == 0.0:  # Check for attachment
                            # once the RNAP has attached, we set the flag has_started to true, so that this condition won't be rechecked.
                            has_started = True
                            start = j
                    else:  # the RNAP has started
                        if self.RNAP_position[j, i] > length:  # The clear sign of detachment of RNAP
                            has_ended = True
                            end = j
                            break
                else:
                    continue
            start_end[i] = [start, end, has_started, has_ended]
        return start_end

    # plotting the trajectory of each RNAP on the axe object provided.
    # each trajectory is ploted individually for cleaness
    def plot_trajectory_RNAP(self, ax):
        final = self.RNAP_position.shape[0]
        # construct time vec
        time = np.arange(0.0, total_time, data_collection_interval)
        start_end = self.find_start_end()

        # plot each line vec
        for j in range(self.RNAP_position.shape[1]):
            start = 0
            end = 0

            if start_end[j, 2]:  # check for start
                start = int(start_end[j, 0])
            else:
                continue

            if start_end[j, 3]:  # check for end
                end = int(start_end[j, 1])
            else:
                end = int(final)
            ax.plot(time[start: end], self.RNAP_position[start:end, j])
        return ax

    # Helper function to get_velocity_not_cleaned()
    # result is not cleaned.
    def get_velocity_not_cleaned(self):
        velocity = np.zeros(self.RNAP_position.shape)

        for i in range(self.RNAP_position.shape[0] - 1):
            for j in range(self.RNAP_position.shape[1]):
                velocity[i, j] = (self.RNAP_position[i, j] - self.RNAP_position[i - 1, j]) / data_collection_interval

        return velocity

    # As the name say, this input the axes object from the Matplotlib and plot the velocity graph on that graph. return the axes
    # plotting the velocity of each RNAP on the axe object provided.
    # velocity graph of each RNAP is plotted separately
    def plot_velocity_RNAP(self, axe, start=0, end=-1):

        if end == -1:
            end = self.RNAP_position.shape[1]

        veloctity_uc = self.get_velocity_not_cleaned()
        time = np.arange(0.0, total_time, data_collection_interval)
        start_end = self.find_start_end()
        final = self.RNAP_position.shape[0]

        # plot each line vec
        for j in range(start, end):
            start = 0
            end = 0

            if start_end[j, 2]:  # check for start
                start = int(start_end[j, 0])
            else:
                continue

            if start_end[j, 3]:  # check for end
                end = int(start_end[j, 1])
            else:
                end = int(final)
            axe.plot(time[start: end], veloctity_uc[start:end, j])
        return axe

    def plot_protein_production_per_second(self, axe):
        protein_ps = np.zeros(int(total_time))
        tot_index = self.protein_production.shape[0]
        n = int(1 / dt)
        for i in range(int(tot_index / n)):
            ps = 0
            for j in range(n):
                ps += self.protein_production[i * n + j]
            protein_ps[i] = ps
        # print all the stuff.

        axe.plot(range(int(total_time)), protein_ps)
        return axe

    def plot_5_and_3_end_amount(self, axe):
        time_list = np.arange(0.0, total_time, data_collection_interval)
        time_list = time_list / 60
        present_5 = self.rnap_loaded - self.rnap_degrading
        present_3 = self.rnap_detached - self.rnap_degraded

        # plot present
        axe.plot(time_list, present_5, 'r--', label='Z5 Present at t')
        axe.plot(time_list, present_3, 'b--', label='Z3 Present at t')

        # plot at
        axe.plot(time_list, self.rnap_loaded, 'r-', label='Z5 Made until t')
        axe.plot(time_list, self.rnap_detached, 'b-', label='Z3 Made until t')
        return axe

    def get_5_and_3_end_amount(self):
        result = np.zeros([4, len(self.rnap_loaded)])
        result[0] = self.rnap_loaded - self.rnap_degrading
        result[1] = self.rnap_detached - self.rnap_degraded
        result[2] = self.rnap_loaded
        result[3] = self.rnap_detached
        return result

# ||---------------------------------------------------------------------||
