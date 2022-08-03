import Model_Protein_Production_OOP as mod
from Model_Protein_Production_OOP import Environment
from Model_Protein_Production_OOP import printProgressBar as printProgressBar
import numpy as np
import matplotlib.pyplot as plt

def uniform_run_repressed(a):
    env = Environment(1/a, include_busty_promoter=False)
    env.DNA.loading_list=mod.uniform_loading_list(mod.total_time, a, 100, 0)
    env.DNA.initiation()
    env.start()
    return env


if __name__ == '__main__':
    # Setting
    mod.testing = True
    mod.group_one = True
    mod.group_two = True
    mod.group_three = True
    mod.group_four = True
    mod.total_time = 500#s 
    alphas = [0.006, 0.033, 0.127]

    # Run
    envs = [None, None, None]
    envs[0] = uniform_run_repressed(alphas[0])
    envs[1] = uniform_run_repressed(alphas[1])
    envs[2] = uniform_run_repressed(alphas[2])

    # File Creation
    path = mod.make_directory_now('bin', 'uniform_runs_repressor')

    # Graphing
    for i in range(3):
        # Path
        name_prom = f'promoter_state_{alphas[i]:.4f}_per_second'.replace('.', '_')

        # Promoter
        fig, ax = plt.subplots(figsize = (10, 5))
        ax.plot(np.arange(0.0, mod.total_time, mod.dt), envs[i].DNA.record_promoter_state)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Promoter State')
        ax.set_title('promoter state graph')
        ax.grid(True)
        ax.set_xlim([0, 250])
        ax.set_ylim([0, 2])
        fig.patch.set_facecolor('xkcd:white')
        fig.savefig(f'{path}/{name_prom}.png')

        target_path = path+f'/loading_list_{alphas[i]:.4f}_per_second'.replace('.', '_')+'.txt'
        np.savetxt(target_path, envs[i].DNA_loading_list, delimiter =', ') 

