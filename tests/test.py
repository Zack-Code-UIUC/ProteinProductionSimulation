import numpy as np
from ProteinProductionSim.variables import k_elong, kRiboLoading, dt, total_time, scaling, length
from ProteinProductionSim.helper.loading_list import LoadingList
from ProteinProductionSim.controller.dna_sim_controller import DNASimController
import matplotlib.pyplot as plt

controller = DNASimController(0.033, promoter_shut_off_time=90, rnap_loading_pattern="uniform")
#controller.total_time = 4
controller.start()
print(controller.env.dna.RNAP_LIST[0].position)
print(controller.env.dna.T_stop)
print(controller.env.dna.detached)
print(controller.env.total_prot)
print(controller.env.dna.loading_list.get())
fig, ax = plt.subplots(figsize=(10, 5))
fig.patch.set_facecolor('xkcd:white')
controller.position.plot(ax)
fig.savefig('position.png')
