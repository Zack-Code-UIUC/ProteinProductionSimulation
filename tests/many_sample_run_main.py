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
        result = [env.DNA.RNAP_LIST[0].attached, env.DNA.RNAP_LIST[0].detach_time, env.get_5_and_3_end_amount()]
    else:
        result = [None]
    return result

def stochastic_run_not_repressed(a):
    env = Environment(1/a, include_busty_promoter=False)
    env.DNA.initiation()
    env.DNA.T_stop = mod.total_time
    env.start()
    if env.DNA.loaded != 0:
        result = [env.DNA.RNAP_LIST[0].attached, env.DNA.RNAP_LIST[0].detach_time, env.get_5_and_3_end_amount()]
    else:
        result = [None]
    return result

if __name__ == "__main__":
    #-----------------=-Setting-----------------
    sample_size = 400
    loading_rates = [0.006, 0.033, 0.127]
    loading_point_low = 0.02
    loading_point_mid = 0.127
    loading_point_high = 1
    experiment_number_stage_1 = 30
    experiment_number_stage_2 = 20
    worker_number = 6

    #--------------------Setup-------------------
    size = experiment_number_stage_1+experiment_number_stage_2+3-1
    data_length = int(mod.total_time/mod.data_collection_interval)
    alphas = np.zeros(size)

    alphas[0] = loading_rates[0]
    alphas[1] = loading_rates[1]
    alphas[2] = loading_rates[2]
    stage_1 = np.linspace(loading_point_low,loading_point_mid, experiment_number_stage_1)
    stage_2 = np.linspace(loading_point_mid,loading_point_high, experiment_number_stage_2)
    alphas[3:experiment_number_stage_1+3]  = stage_1
    alphas[experiment_number_stage_1+3:experiment_number_stage_2+experiment_number_stage_1+3] = stage_2[1:]
    finished = np.zeros([2,size])
    finished_time_sum = np.zeros([2,size])
    five_three_sum = np.zeros([2, size, 4, data_length])
    pool = multiprocessing.Pool(worker_number)

    #-----Stochastic Run With Repression--------

    t_i = time.perf_counter()

    printProgressBar(0, sample_size*2)
    for i in range(sample_size):
        
        result = pool.map(stochastic_run_repressed, alphas)
        
        for j in range(size):
            if result[j][0] == False:
                finished[0,j] +=1
                finished_time_sum[0,j] += result[j][1]
                five_three_sum[0, j] = five_three_sum[0, j]+np.array(result[j][2])
            elif result[j][0] == True:
                five_three_sum[0, j] = five_three_sum[0, j]+np.array(result[j][2])

        printProgressBar(i, sample_size*2)
    
    #-----Stochastic Run Without Repression--------

    for i in range(sample_size):

        result = pool.map(stochastic_run_not_repressed, alphas)
        
        for j in range(size):
            if result[j][0] == False:
                finished[1,j] +=1
                finished_time_sum[1,j] += result[j][1]
                five_three_sum[1, j] = five_three_sum[1, j]+np.array(result[j][2])
            elif result[j][0] == True:
                five_three_sum[1, j] = five_three_sum[1, j]+np.array(result[j][2])

        printProgressBar(i+sample_size, sample_size*2)

    printProgressBar(sample_size*2, sample_size*2)
    t_f = time.perf_counter()
    print(f'Simulation Runtime: {(t_f-t_i)/60} min.')

    #--------------Data Analysis-----------------
    t_i = time.perf_counter()
        
    finished_rate_repressed = finished[0]/sample_size
    finished_time_average_repressed =np.zeros(size)
    for i in range(size):
        if finished[0,i] == 0:
            finished_time_average_repressed[i] = 0
        else:
            finished_time_average_repressed[i] = finished_time_sum[0,i]/finished[0,i]
    five_three_average_repressed = five_three_sum[0]/sample_size

    finished_rate_not_repressed = finished[1]/sample_size
    finished_time_average_not_repressed =np.zeros(size)
    for i in range(size):
        if finished[1,i] == 0:
            finished_time_average_not_repressed[i] = 0
        else:
            finished_time_average_not_repressed[i] = finished_time_sum[1,i]/finished[1,i]
    five_three_average_not_repressed = five_three_sum[1]/sample_size

    t_f = time.perf_counter()
    print(f'Analysis Runtime: {t_f-t_i} s.')
    #------File Creation and Storing Setting---------
    t_i = time.perf_counter()
    # create a folder to store the data
    path = mod.make_directory_now('bin', 'many_sample_runs')

    # write run setting:
    setting_path = path+'\\setting.txt'
    setting = open(setting_path, 'w')
    setting.write(f'total_time: {mod.total_time}\n')
    setting.write(f'sample_size: {sample_size}\n')
    setting.write(f'worker_number: {worker_number}\n')
    setting.write(f'alphas: {str(alphas)}\n')

    # write all array data
    target_path = path+'\\finished_rate_repressed.txt'
    np.savetxt(target_path, finished_rate_repressed, delimiter =', ') 
    target_path = path+'\\finished_time_average_repressed.txt'
    np.savetxt(target_path, finished_time_average_repressed, delimiter =', ') 
    target_path = path+'\\five_three_average_repressed.npy'
    np.save(target_path, five_three_average_repressed)
     
    target_path = path+'\\finished_rate_not_repressed.txt'
    np.savetxt(target_path, finished_rate_not_repressed, delimiter =', ') 
    target_path = path+'\\finished_time_average_not_repressed.txt'
    np.savetxt(target_path, finished_time_average_not_repressed, delimiter =', ') 
    target_path = path+'\\five_three_average_not_repressed.npy'
    np.save(target_path, five_three_average_not_repressed) 

    t_f = time.perf_counter()
    print(f'File Writing Runtime: {t_f-t_i} s.')
    #-----------Data Visualization-----------------

    t_i = time.perf_counter()
    #----Repressed----
    # finish rate graph
    fig, ax = plt.subplots(figsize = (12, 5))
    x=np.linspace(-0.1,loading_point_high+0.2,100)
    y=np.zeros(100)+1
    ax.plot(x,y,'--k')
    ax.set_xlabel(r'$\alpha$ $\mathregular{(s^{-1})}$', fontsize=15)
    ax.set_ylabel('Percentage of RNAP finished (bp/s)', fontsize=13)
    ax.plot(alphas[3:],finished_rate_repressed[3:],color='blue', marker='^', linestyle='solid',linewidth=2, markersize=8,label='simulation')
    ax.set_xticks(ticks = np.arange(0, 1.2, step = 0.05))
    ax.set_yticks(ticks = np.arange(0, 1, step=0.1))
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.2])
    ax.set_title('The Finish Rate of the First RNAP With Repressor On at 90s')
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig(f'{path}\\finish_rate_repressed.png')

    # finish time graph
    fig, ax = plt.subplots(figsize = (12, 5))
    ax.set_xlabel(r'$\alpha$ $\mathregular{(s^{-1})}$', fontsize=15)
    ax.set_ylabel('Average Finish-Time(s)', fontsize=15)
    ax.set_title('Average Finish-Time of the First RNAP With Repressor On at 90s')
    ax.plot(alphas[3:],finished_time_average_repressed[3:],color='blue', marker='^', linestyle='solid',linewidth=2, markersize=6,label='simulation')
    ax.set_xticks(ticks = np.arange(0, 1.2, step = 0.05))
    ax.set_xlim([0, 1])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig(f'{path}\\average_finish_time_repressed.png')

    # five three graph of the low, intermediate and high loading rate. 
    x=np.arange(0.0, mod.total_time, mod.data_collection_interval)/60
    for i in range(3):
        fig, axe = plt.subplots(figsize = (12, 8))
        fig.patch.set_facecolor('xkcd:white')
        y = five_three_average_repressed[i]
        axe.plot(x, y[0], 'r--', label = 'Z5 Present at t')
        axe.plot(x, y[1], 'b--', label = 'Z3 Present at t')

        # plot at
        axe.plot(x, y[2], 'r-', label = 'Z5 Made until t')
        axe.plot(x, y[3], 'b-', label = 'Z3 Made until t')
        axe.set_xlabel('Time[min]', fontsize = 15)
        axe.set_ylabel('Number of RNA Ends', fontsize = 15)
        axe.set_title('5\' and 3\' End amount Versus Time With Repressor On at 90s')
        axe.legend()
        axe.grid(True)
        name = f'fiveandthree_repressed_{alphas[i]:.4f}_per_second'.replace('.', '_')
        fig.savefig(f'{path}\\{name}.png')

    #----Not Repressed----
    # finish rate graph
    fig, ax = plt.subplots(figsize = (12, 5))
    x=np.linspace(-0.1,loading_point_high+0.2,100)
    y=np.zeros(100)+1
    ax.plot(x,y,'--k')
    ax.set_xlabel(r'$\alpha$ $\mathregular{(s^{-1})}$', fontsize=15)
    ax.set_ylabel('Percentage of RNAP finished (bp/s)', fontsize=13)
    ax.plot(alphas[3:],finished_rate_not_repressed[3:],color='blue', marker='^', linestyle='solid',linewidth=2, markersize=8,label='simulation')
    ax.set_xticks(ticks = np.arange(0, 1.2, step = 0.05))
    ax.set_yticks(ticks = np.arange(0, 1, step=0.1))
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.2])
    ax.set_title('The Finish Rate of the First RNAP With No Repression')
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig(f'{path}\\finish_rate_not_repressed.png')

    # finish time graph
    fig, ax = plt.subplots(figsize = (12, 5))
    ax.set_xlabel(r'$\alpha$ $\mathregular{(s^{-1})}$', fontsize=15)
    ax.set_ylabel('Average Finish-Time(s)', fontsize=15)
    ax.set_title('Average Finish-Time of the First RNAP With No Repression')
    ax.plot(alphas[3:],finished_time_average_not_repressed[3:],color='blue', marker='^', linestyle='solid',linewidth=2, markersize=6,label='simulation')
    ax.set_xticks(ticks = np.arange(0, 1.2, step = 0.05))
    ax.set_xlim([0, 1])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig(f'{path}\\average_finish_time_not_repressed.png')

    # five three graph of the low, intermediate and high loading rate. 
    x=np.arange(0.0, mod.total_time, mod.data_collection_interval)/60
    for i in range(3):
        fig, axe = plt.subplots(figsize = (12, 8))
        fig.patch.set_facecolor('xkcd:white')
        y = five_three_average_not_repressed[i]
        axe.plot(x, y[0], 'r--', label = 'Z5 Present at t')
        axe.plot(x, y[1], 'b--', label = 'Z3 Present at t')

        # plot at
        axe.plot(x, y[2], 'r-', label = 'Z5 Made until t')
        axe.plot(x, y[3], 'b-', label = 'Z3 Made until t')
        axe.set_xlabel('Time[min]', fontsize = 15)
        axe.set_ylabel('Number of RNA Ends', fontsize = 15)
        axe.set_title('5\' and 3\' End amount Versus Time With Repressor On at 90s')
        axe.legend()
        axe.grid(True)
        name = f'fiveandthree_not_repressed_{alphas[i]:.4f}_per_second'.replace('.', '_')
        fig.savefig(f'{path}\\{name}.png')

    t_f = time.perf_counter()
    print(f'Graphing Runtime: {t_f-t_i} s.')