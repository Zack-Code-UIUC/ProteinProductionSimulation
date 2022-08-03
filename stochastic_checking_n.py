import Model_Protein_Production_OOP as mod
from Model_Protein_Production_OOP import Environment
from Model_Protein_Production_OOP import printProgressBar as printProgressBar
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import time


def stochastic_run_repressed(a):
    env = Environment(1/a, include_busty_promoter=False)
    env.DNA.initiation()
    env.start()
    if env.DNA.loaded != 0:
        result = [env.DNA.RNAP_LIST[0].attached, env.DNA.RNAP_LIST[0].detach_time, env.get_5_and_3_end_amount(), env.DNA.loaded]
    else:
        result = [None]
    return result

def stochastic_run_not_repressed(a):
    env = Environment(1/a, include_busty_promoter=False)
    env.DNA.T_stop=mod.total_time
    env.DNA.initiation()
    env.start()
    if env.DNA.loaded != 0:
        result = [env.DNA.RNAP_LIST[0].attached, env.DNA.RNAP_LIST[0].detach_time, env.get_5_and_3_end_amount(), env.DNA.loaded]
    else:
        result = [None]
    return result

if __name__ == "__main__":
    #-----------------=-Setting-----------------
    sample_size = 5
    loading_rates = [0.006, 0.033, 0.127]
    loading_point_low = 0.02
    loading_point_high = 0.127
    experiment_number_stage = 50
    worker_number = 6
    #--------------------Setup-------------------
    size = experiment_number_stage+3-1
    data_length = int(mod.total_time/mod.data_collection_interval)
    alphas = np.zeros(size)

    alphas[0] = loading_rates[0]
    alphas[1] = loading_rates[1]
    alphas[2] = loading_rates[2]
    stage = np.linspace(loading_point_low,loading_point_high, experiment_number_stage)
    alphas[3:experiment_number_stage+3]  = stage[:-1]
    #print(alphas)
    finished = np.zeros([size])
    finished_time_sum = np.zeros([size])
    five_three_sum = np.zeros([size, 4, data_length])
    finished_rnap_number = np.zeros([2, sample_size*size])
    pool = multiprocessing.Pool(worker_number)

    #-----Stochastic Run With Repression--------

    t_i = time.perf_counter()
    printProgressBar(0, sample_size)
    for i in range(sample_size):
        result = pool.map(stochastic_run_repressed, alphas)
        for j in range(size):
            index = i*sample_size+j
            if result[j][0] == False:
                finished[j] +=1
                finished_time_sum[j] += result[j][1]
                five_three_sum[j] = five_three_sum[j]+np.array(result[j][2])
                finished_rnap_number[0, index] = 1
                finished_rnap_number[1, index] = result[j][3]
            elif result[j][0] == True:
                five_three_sum[j] = five_three_sum[j]+np.array(result[j][2])
                finished_rnap_number[0, index] = 0
                finished_rnap_number[1, index] = result[j][3]

        printProgressBar(i+1, sample_size)
    t_f = time.perf_counter()
    print(f'Simulation Runtime: {(t_f-t_i)/60} min.')
    #--------------Data Analysis-----------------
    t_i = time.perf_counter()
        
    finished_rate_repressed = finished/sample_size
    finished_time_average_repressed =np.zeros(size)
    for i in range(size):
        if finished[i] == 0:
            finished_time_average_repressed[i] = 0
        else:
            finished_time_average_repressed[i] = finished_time_sum[i]/finished[i]
    five_three_average_repressed = five_three_sum/sample_size

    tiles_number = 50 + 1
    tiles = range(tiles_number)
    amount_each_tile = np.zeros(tiles_number)
    finished_each_tile = np.zeros(tiles_number)

    for i in range(sample_size*size):
        rnap_amount = int(finished_rnap_number[1, i])
        amount_each_tile[rnap_amount] +=1
        if finished_rnap_number[0, i] == 1:
            finished_each_tile[rnap_amount] +=1
    finish_rate_data = []
    finish_rate_n = []
    for i in range(tiles_number):
        if amount_each_tile[i] != 0 and i != 0:
            finish_rate_data.append(finished_each_tile[i]/amount_each_tile[i])
            finish_rate_n.append(i)

    t_f = time.perf_counter()
    print(f'Analysis Runtime: {t_f-t_i} s.')

    #------File Creation and Storing Setting---------
    t_i = time.perf_counter()
    # create a folder to store the data
    path = mod.make_directory_now('bin', 'stochastic_checking_n')

    # write run setting:
    setting_path = path+'/setting.txt'
    setting = open(setting_path, 'w')
    setting.write(f'total_time: {mod.total_time}\n')
    setting.write(f'sample_size: {sample_size}\n')
    setting.write(f'worker_number: {worker_number}\n')
    setting.write(f'alphas: {str(alphas)}\n')

    # write all array data
    target_path = path+'/finished_rate_repressed.txt'
    np.savetxt(target_path, finished_rate_repressed, delimiter =', ') 
    target_path = path+'/finished_time_average_repressed.txt'
    np.savetxt(target_path, finished_time_average_repressed, delimiter =', ') 
    target_path = path+'/five_three_average_repressed.npy'
    np.save(target_path, five_three_average_repressed)
    target_path = path+'/finished_rnap_number.npy'
    np.save(target_path, finished_rnap_number) 

    t_f = time.perf_counter()
    print(f'File Writing Runtime: {t_f-t_i} s.')

    #-----------Data Visualization-----------------

    t_i = time.perf_counter()

    # finish rate graph
    fig, ax = plt.subplots(figsize = (12, 5))
    ax.set_xlabel('Number of RNAP loaded', fontsize=15)
    ax.set_ylabel('Percentage of Cases which the First RNAP finished', fontsize=13)
    ax.plot(finish_rate_n,finish_rate_data,color='blue', marker='^', linestyle='None',linewidth=2, markersize=8,label='simulation')
    ax.set_title('The Finish Rate of the First RNAP versus RNAPs loaded With Repressor On at 90s')
    fig.patch.set_facecolor('xkcd:white')
    ax.set_xticks(np.arange(np.min(finish_rate_n), np.max(finish_rate_n)+1, 1.0))
    fig.savefig(f'{path}/finish_rate_versus_n.png')

    t_f = time.perf_counter()
    print(f'Graphing Runtime: {t_f-t_i} s.')