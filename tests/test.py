import numpy as np
from ProteinProductionSim.variables import k_elong, kRiboLoading, dt, total_time, scaling
from ProteinProductionSim.helper.loading_list import LoadingList

l = LoadingList(None, 1000, 1/40, if_stochastic= True, if_bursty=False)
print(l.get())
print(l.get_average_loading_interval())
l.increment()
print(l.if_can_load(30))
print(l.get_location())
l.check_for_repeat()
l.trim(500)
print(l.get())
print(scaling(total_time))
