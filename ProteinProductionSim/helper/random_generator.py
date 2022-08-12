"""
===================
random_generator.py
===================

This is the file for all the random generator.
"""


from numpy import random as rand
import math


def exponential_generator(rate):
    """
    This method generate a float value based on a exponential distribution.

    Parameters
    ----------
    rate : float
        this is the rate of the exponential distribution or 1/mean of the pdf.

    Returns
    -------
    float
    """
    return rand.exponential(scale=1/rate)


def stepwise_exponential_generator(m1, m2, t_crit):
    """
    This method generate a float value based on a step-wise random distribution.

    The pdf is that of two combined exponential pdf.

    Parameters
    ----------
    m1 : float
        this represents the mean of the front exponential distribution
    m2 : float
        this represents the mean of the back exponential distribution
    t_crit : float
        this is the cutoff value between the two exponential distribution

    Returns
    ------
    float
        the end result
    """
    portion1 = 1-math.exp(-1*t_crit/m1)
    portion2 = 1-portion1
    choice = rand.choice([1,2], p = [portion1, portion2])
    if choice == 1:
        passed = False
        while not passed:
            result =rand.exponential(scale = m1)
            if result <= t_crit:
                passed = True
    else:
        passed = False
        while not passed:
            result =rand.exponential(scale = m2)
            if result >= t_crit:
                passed = True
    return result


def binary_generator(probability):
    """
    This method generate a true or false value randomly based on the provided probability

    Parameters
    ----------
    probability : float
        the probability for the result to be true

    Returns
    -------
    bool
        True of False
    """
    return rand.choice(a=[True, False], p=[1 - probability, probability])