import numpy as np
if True:
    loading_rates = [0.006, 0.033, 0.127]
    loading_point_low = 0.02
    loading_point_high = 0.127
    experiment_number_stage = 50
    worker_number = 6
    #--------------------Setup-------------------
    size = experiment_number_stage+3-1
    alphas = np.zeros(size)

    alphas[0] = loading_rates[0]
    alphas[1] = loading_rates[1]
    alphas[2] = loading_rates[2]
    stage = np.linspace(loading_point_low,loading_point_high, experiment_number_stage)
    alphas[3:experiment_number_stage+3]  = stage[:-1]
    print(alphas)