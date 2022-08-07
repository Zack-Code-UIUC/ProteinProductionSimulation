# Helper Functions
# ||---------------------------------------------------------------------||
# Supercoiling
def S(R):
    return (gamma * R)


def n_dependence_cubic_3(x):
    # return 1+0.778753*(x-1)+3.3249*(x-1)**2+0.379478*(x-3)**3
    return 4.17691 - 6.16402 * x + 2.73771 * x ** 2 + 0.249398 * x ** 3


def phi(RNAP_list, detached_rnap_amount):
    size = len(RNAP_list) - detached_rnap_amount
    if size == 0:
        return None
    PHI = np.zeros(size + 1)

    for i in range(size + 1):
        j = i + detached_rnap_amount
        if i == 0:
            PHI[i] = S(RNAP_list[j].position)
            continue
        elif i == size:
            PHI[i] = 0  # S(RNAP_list[j-1].position)
            continue
        else:
            PHI[i] = S(RNAP_list[j - 1].position - RNAP_list[j].position)

    return PHI


def torque(phi):
    size = np.size(phi) - 1
    torq = np.zeros(size)
    n = n_dependence_cubic_3(size)
    for i in range(size):
        torq[i] = -tau_0 * n * (phi[i] - phi[i + 1])

    return torq


def velocity(tau, detached_rnap_amount):
    size = np.size(tau) + detached_rnap_amount
    velo = np.zeros(size)

    for i in range(detached_rnap_amount):
        velo[i] = 0

    for i in range(np.size(tau)):
        j = i + detached_rnap_amount
        if tau[i] > 1.5 * tau_c:
            velo[j] = 0
        elif tau[i] < -1.5 * tau_c:
            velo[j] = 2 * v_0
        else:
            velo[j] = 2 * v_0 / (1 + np.exp(2 * (tau[i] / tau_c) ** 3))
    return velo


def supercoil_velocity(RNAP_list):
    detached_rnap_amount = 0

    for rnap in RNAP_list:
        if not rnap.attached:
            detached_rnap_amount += 1

    if detached_rnap_amount == len(RNAP_list):
        return np.zeros(len(RNAP_list))

    PHI = phi(RNAP_list, detached_rnap_amount)
    torq = torque(PHI)
    velo = velocity(torq, detached_rnap_amount)
    return velo


# ||---------------------------------------------------------------------||
# RNAP_loading_list

# RNAP_time_list(arr, load_rate) is used to generate loading attempt time.
def RNAP_loading_list(arr, load_rate):
    t_loading = []
    time_last = 0
    for pair in arr:
        if pair[0]:
            loading_list = load_list(pair[1] - time_last, load_rate)
            acc_t = 0
            for t in loading_list:
                t_loading.append(t + time_last + acc_t)
                acc_t += t
        time_last = pair[1]
    # print(len(t_loading))
    # additional normalization is performed to clean any loading attempt that is too close to the other set temporally.
    # so that within one dt, only one attempt is performed. Such attempt would obviously fail.
    elem = 0
    to_remove = []
    for count, t in enumerate(t_loading):
        if count != 0 and not (t in to_remove):
            if t - elem < dt:
                to_remove.append(t)
            else:
                elem = t
    for item in to_remove:
        t_loading.remove(item)

    return t_loading


def RIBO_loading_list(duration, load_rate):
    loading_list = load_list(duration, load_rate)
    t_loading = []
    acc_t = 0
    for t in loading_list:
        t_loading.append(t + acc_t)
        acc_t += t
    return t_loading


def load_list(duration, load_rate):
    ave = 1 / load_rate
    t = 0
    t_slots = []

    while t <= duration:
        add_time = rand.exponential(scale=ave)
        if add_time < duration - t:
            t_slots.append(add_time)
            t += add_time
        else:
            break

    return t_slots


# testing
# arr = promoter_time_list(7000,5,3)
# print(arr)
# print(RNAP_loading_list(arr, 1/15))
# print(load_list(20,1/10))
# print(sum(load_list(10,1/10)))
# print(RIBO_loading_list(30, 1/15))

# ||---------------------------------------------------------------------||
# promoter_time_list
# This function is used for class DNAstrand initiation. It chucks out a 1-D np array with lists as its elements.
# Each list represents a time slot with the status of the promoter and the end time of the status.
# For instance, [True, 13.00], means that the status of the promoter is On, and the status will remain On until 13.00 time_point

# The time is generated following exponential waiting time distribution.

def promoter_time_list(time_tot, on_time, off_time):
    t_slots = promoter_time_list_non_cumulative(time_tot, on_time, off_time)
    tot_t = 0
    for slot in t_slots:
        slot[1] = tot_t + slot[1]
        tot_t = slot[1]
    return t_slots


def promoter_time_list_non_cumulative(time_tot, on_time, off_time):
    t = 0
    # True indicates a On Promoter. And False indicates an Off Promoter
    i = False  # This value here, represent the initial state of the Promoter
    t_slots = []
    while t <= time_tot:
        if i:
            add_time = rand.exponential(scale=off_time)
            if add_time < time_tot - t:
                i = not i
                t_slots.append([False, add_time])
                t += add_time
            else:
                t_slots.append([False, time_tot - t])
                break
        else:
            add_time = rand.exponential(scale=on_time)
            if add_time < time_tot - t:
                i = not i
                t_slots.append([True, add_time])
                t += add_time
            else:
                t_slots.append([True, time_tot - t])
                break

    return t_slots


# Testing
# arr = promoter_time_list(100,5,3)
# print(arr)

def average_ON_OFF_Time(slots):
    total_ON_t = 0
    total_OFF_t = 0
    on_n = 0
    off_n = 0
    for slot in slots:
        if slot[0]:
            on_n += 1
            total_ON_t = total_ON_t + slot[1]
        else:
            off_n += 1
            total_OFF_t = total_OFF_t + slot[1]
    return [total_ON_t / on_n, total_OFF_t / off_n, on_n, off_n]


# arr2 = promoter_time_list_non_cumulative(100,5,3)
# print(average_ON_OFF_Time(arr2))
# ||---------------------------------------------------------------------||
# Additional Promoter Load Lists

# For this loading list generates a uniform non-busty
def uniform_promoter_time_list(time_tot, on_time, off_time):
    t = 0
    # True indicates a On Promoter. And False indicates an Off Promoter
    i = False  # This value here, represent the initial state of the Promoter
    t_slots = []
    while t <= time_tot:
        if i:
            add_time = off_time
            if add_time < time_tot - t:
                i = not i
                t_slots.append([False, add_time])
                t += add_time
            else:
                t_slots.append([False, time_tot - t])
                break
        else:
            add_time = on_time
            if add_time < time_tot - t:
                i = not i
                t_slots.append([True, add_time])
                t += add_time
            else:
                t_slots.append([True, time_tot - t])
                break
    tot_t = 0
    for slot in t_slots:
        slot[1] = tot_t + slot[1]
        tot_t = slot[1]
    return t_slots


def uniform_loading_list(tot_time, loading_rate, repeating_loading, off_interval):
    loading_time = 1 / loading_rate

    load_list = [0]
    t = 0.0
    prom_flag = True
    have_loaded = 1
    while t <= tot_time:
        if prom_flag:
            if t + loading_time >= tot_time:
                break
            if have_loaded < repeating_loading - 1:
                load_list.append(t + loading_time)
                t += loading_time
                have_loaded += 1
            elif have_loaded == repeating_loading - 1:
                load_list.append(loading_time + t)
                t += loading_time
                have_loaded = 0
                prom_flag = False
        else:
            t += off_interval
            prom_flag = True

    return load_list


# ||---------------------------------------------------------------------||
# Custom Random Generator
# this is used for generator two phase step-wise random distribution

def stepwise_exponential_generator(m1, m2, tcrit):
    # factor = math.exp(-1*tcrit/m1)
    portion1 = 1 - math.exp(-1 * tcrit / m1)
    portion2 = 1 - portion1
    choice = rand.choice([1, 2], p=[portion1, portion2])
    if choice == 1:
        passed = False
        while not passed:
            result = rand.exponential(scale=m1)
            if result <= tcrit:
                passed = True
    else:
        passed = False
        while not passed:
            result = rand.exponential(scale=m2)
            if result >= tcrit:
                passed = True
    return result


# Variable for this function: 120, 90, 102.
# fig, ax = plt.subplots(1, sharey=True, tight_layout=True)

# stepwise_exponential_generator(20, 10, 5)
# sampling = np.zeros(10000)
# for i in range(10000):
# sampling[i] = stepwise_exponential_generator(90,45, 102)

# ax.hist(sampling, bins = 1000)
# ax.set_title('Histogram of the Stepwise Exponential Random Generator: 10,000 Attempt')
# ax.grid(True);
# fig.savefig('Histogram_Stepwise_Exp_Rand_Generator_10,000_Attempt')
# ||---------------------------------------------------------------------||
# Progress Bar

def printProgressBar_RNAP_prot(iteration, total, total_RNAP, total_prot):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        total_RNAP  - Required  : total number of RNAPs that have been register in the DNA instance (Int)
        total_prot  - Required  : total number of proteins that have been produced (Int)
    """
    prefix = 'Progress:'
    suffix = 'Completed'
    UP = "\x1B[2A"
    CLR = "\x1B[0K"
    fill = '█'
    length = 50
    decimals = 2

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{CLR}{prefix} |{bar}| {percent}% {suffix}',
          '\n{:<15} {:<8} {:<15} {:<8}\n'.format('RNAP Number =', str(total_RNAP), 'Protein Number =', str(total_prot)),
          end=UP)

    if iteration == total:
        print()
        print()


def printProgressBar(iteration, total):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
    """
    prefix = 'Progress:'
    suffix = 'Completed'
    fill = '█'
    length = 50
    decimals = 2

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'{prefix} |{bar}| {percent}% {suffix}', end='\r')

    if iteration == total:
        print()


# ||---------------------------------------------------------------------||
# Make Directory based on time and date
def make_directory_now(direc, name):
    '''
    Create a directory based on the target directory and inputed name. A date and time string will be attached to the back of the name.

    Parameters:
    direc(string): target directory
    name(string) : directory name that need to be created.

    Return:
    string: the path of the create directory
    '''
    # get current working directory
    cwd = os.getcwd()
    # construct the date and time string
    now = datetime.now()
    date_string = now.strftime("%Y_%m_%d_%H_%M")
    # make the directory name
    directory_name = name + '_' + date_string
    # construct the whole path
    path_directory = cwd + '/' + direc + '/' + directory_name
    try:
        makedirs(path_directory)
    except FileExistsError:
        pass
    return path_directory
# ||---------------------------------------------------------------------||
