"""
===============
loading_list.py
===============

This helper file contains all the helper function related to the loading list for both RNAPs and Ribosomes.

Loading lists are expressed in terms of python lists. Each element represents a time point(in integer index). The
order of the element is always ascending for each comparison. All the loading list are subclasses of the
DataContainer.

A loading list can be either stochastic or uniform.

A loading list can be used for either RNAP(DNA) or Ribosomes(mRNA). For RNAP loading, the loading list can also be
bursty or non-bursty, which means that it also have a promoter list aside from the list. It must be noted that his
promoter list is different from the actual "real" promoter list, as that will be record by the controller using
other DataContainer.
"""

from interface import DataContainer
from .random_generator import exponential_generator

# Helper functions for the helper class


def _stochastic_cumulative_array_generator(duration, rate):
    """
    This method generate a stochastic loading array. The elements are accumulative.
    Parameters
    ----------
    duration : int
    rate : float

    Returns
    -------
    list[float]
    """
    t = 0.0
    t_slots = [0]
    while t <= duration:
        add_time = exponential_generator(rate)
        if add_time < float(duration) - t:
            t += add_time
            t_slots.append(t)
        else:
            break
    return t_slots


def _stochastic_noncumulative_array_generator(duration, rate):
    """
    This method generate a stochastic loading array. The elements are noncumulative.
    Parameters
    ----------
    duration : int
    rate : float

    Returns
    -------
    list[float]
    """
    t = 0.0
    t_slots = [0]
    while t <= duration:
        add_time = exponential_generator(rate)
        if add_time < float(duration) - t:
            t_slots.append(t)
            t += add_time
        else:
            break
    return t_slots


def _uniform_cumulative_array_generator(duration, rate):
    """
    This method generate a uniform loading array. The elements are noncumulative.
    Parameters
    ----------
    duration : int
    rate : float

    Returns
    -------
    list[float]
    """
    interval = 1/rate
    t = 0.0
    t_slots = [0]
    while t <= duration:
        add_time = interval
        if add_time < float(duration) - t:
            t += add_time
            t_slots.append(t)
        else:
            break
    return t_slots


def _uniform_noncumulative_array_generator(duration, rate):
    """
    This method generate a loading array. The elements are noncumulative.
    Parameters
    ----------
    duration : int
    rate : float

    Returns
    -------
    list[float]
    """
    interval = 1/rate
    t = 0.0
    t_slots = [0]
    while t <= duration:
        add_time = interval
        if add_time < float(duration) - t:
            t_slots.append(t)
            t += add_time
        else:
            break
    return t_slots


def _promoter_array_generator(duration, tau_on, tau_off):
    duration = float(duration)
    t = 0
    i = False
    t_slots = []
    while t <= duration:
        if i:
            add_time = exponential_generator(1/tau_off)
            if add_time < duration - t:
                i = not i
                t_slots.append([False, add_time])
                t += add_time
            else:
                t_slots.append([False, duration - t])
                break
        else:
            add_time = exponential_generator(1/tau_on)
            if add_time < duration - t:
                i = not i
                t_slots.append([True, add_time])
                t += add_time
            else:
                t_slots.append([True, duration - t])
                break

    return t_slots


def _stochastic_bursty_array_generator(duration, rate, tau_off=143.0, tau_loading=2.2):
    """
    This method generate a loading list that is both bursty and stochastic.
    """
    tau_on = rate * tau_loading * tau_off / (1.0 - rate * tau_loading)
    promoter_list = _promoter_array_generator(duration, tau_on, tau_off)
    t_loading = []
    for pair in promoter_list:
        if pair[0]:
            loading_list = _stochastic_noncumulative_array_generator(pair[1], 1/tau_loading)
            for t in loading_list:
                t_loading.append(t)

    for i in range(len(t_loading)-1):
        t_loading[i+1] = t_loading[i] + t_loading[i]
    return t_loading, promoter_list


def _uniform_bursty_array_generator(duration, rate, tau_off=143.0, tau_loading=2.2):
    """
    This method generate a loading list that is both bursty and stochastic.
    """
    tau_on = rate * tau_loading * tau_off / (1.0 - rate * tau_loading)
    promoter_list = _promoter_array_generator(duration, tau_on, tau_off)
    t_loading = []
    for pair in promoter_list:
        if pair[0]:
            loading_list = _uniform_noncumulative_array_generator(pair[1], 1 / tau_loading)
            for t in loading_list:
                t_loading.append(t)

    for i in range(len(t_loading) - 1):
        t_loading[i + 1] = t_loading[i] + t_loading[i]
    return t_loading, promoter_list


class LoadingList(DataContainer):
    """
    This is a container class for the loading list.

    Attributes
    ----------
    location : int
        this represents the current index of the loading list.
    arr : numpy array of int
        this represents the loading list, each element is in index form.
    """
    def __init__(self, parent, duration, rate, if_stochastic=False, if_bursty=False):
        super().__init__(parent)
        self.location = 0
        self._arr = None
        self._promoter_list = None
        match if_bursty:
            case False:
                match if_stochastic:
                    case True:
                        self._arr = _stochastic_cumulative_array_generator(duration, rate)
                    case False:
                        self._arr = _uniform_cumulative_array_generator(duration, rate)
            case True:
                self._arr, self.promoter_list = _stochastic_bursty_array_generator(duration, rate)

        self.length = len(self._arr)
        self.arr = [int(self._arr[i]) for i in range(self.length)]
        self.length = len(self.arr)
        self._remove_duplicate()
        self.dumped = False

    def init(self):
        pass

    def check_for_repeat(self):
        for i in range(len(self.arr)-1):
            if self.arr[i] == self.arr[i+1]:
                print("Repeat!")

    def _remove_duplicate(self):
        self.arr = [i for n, i in enumerate(self.arr) if i not in self.arr[:n]]
        self.length = len(self.arr)

    def log(self, **kwargs):
        pass

    def get(self):
        return self.arr

    def get_location(self):
        return self.location

    def get_current(self):
        return self.arr[self.location]

    def get_promoter(self):
        return self.promoter_list

    def increment(self):
        self.location += 1

    def get_average_loading_interval(self):
        if self.length == 0:
            return 0
        tot = 0
        for i in range(self.length-1):
            tot += self.arr[i + 1] - self.arr[i]
        return tot / (self.length-1)

    def if_empty(self):
        return self.location >= self.length

    def dump(self):
        self.length = 0
        self.arr = []

    def if_can_load(self, time_index, base_time=0):
        if not self.if_empty() and time_index == self.get_current() + base_time:
            self.increment()
            return True
        else:
            return False

    def trim(self, t_stop):
        self.arr = [i for i in self.arr if i <= t_stop]
        self.length = len(self.arr)

    def get_length(self):
        return self.length
