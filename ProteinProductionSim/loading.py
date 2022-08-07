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


