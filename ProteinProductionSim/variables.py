# Parameters

# Environment variable
import numpy as np

total_time = 300 # simulation time second
dt = 1/30 # s
data_collection_interval = 1/10  # second
stage_per_collection = int(data_collection_interval/dt)
data_smoothing_interval = 5  # s

# Basic Coefficients
length = 3072   # Length of lacZ gene or 3075 bps
# proteinLL = 1200 #s

# Busty Promoter Coefficient
t_on = 7.8 #s
stop=np.array([total_time,200,400])
tau_off= 143 # s 1/k_off
tau_loading = 1/0.033 # s 1/k_loading
list_interval_loading = [3.5, 7, 15, 30, 500] # s 1/k_on

# Site-Specific Pausing Coefficient
pauseProfile = ("flat", "OnepauseAbs", "TwopauseAbs")
pauseSite = np.array((1500, 2500))
pauseDuration = np.array((10, 15))  # s
pauseProb = 0.8

# RNAP Variables
RNAP_size = 35 # bps
k_elong = 30.5
mRNALL = 90 # s Mean mRNA lifetime
proteinLL = 1200 #s


# Ribosome Variables
RIBO_size = 30 #nts
ribo_loading_interval = 90
tau_mRNA = 90
kRiboLoading = 0.2 # sec^-1
initiation_nt = 33 # nts
m1=90
m2=45
t_crit=102

# mRNA degradation
list_times_degradation = (7, 15, 30, 60) #s

# Supercoiling Variables
gamma = 0.01 # supercoiling constant
v_0 = 30.5
tau_c = 11
tau_0 = 0.386


# Scaling
multiplier = int(1/dt)


def scaling(variables):
    return int(multiplier * variables)
