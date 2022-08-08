
#||=====================================================================||
# Simulation

# Run Setting Adjustment
if __name__ == '__main__':
    testing = True
    group_one = True
    group_two = True
    group_three = True
    group_four = True
    total_time = 500#s
    #data_collection_interval = dt
    env = Environment(1/0.127, include_busty_promoter=False)
    env.DNA.loading_list = uniform_loading_list(total_time, 0.12832653061224492, 80, 0)
    env.DNA.T_stop = 90
    env.DNA.initiation()
    env.start()

#||=====================================================================||
    #Graphing
    print("Graphing Starts.")
    start_time = time.perf_counter()
    time_list = np.arange(0.0, total_time, data_collection_interval)
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('xkcd:white')
    ax.set_title('RNAP Amount')
    ax.plot(time_list, env.RNAP_amount)
    ax.grid(True)
    fig.savefig('rnap_amount.png')

#||---------------------------------------------------------------------||

    fig, ax = plt.subplots(figsize = (10, 5))
    fig.patch.set_facecolor('xkcd:white')
    env.plot_trajectory_RNAP(ax)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Position [bps]')
    ax.set_title('RNAP Position Plot')
    ax.set_xlim([-10, 240])
    #ax.axis([-5,233.3, -100, length+200])
    ax.grid(True)
    fig.savefig('position.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    time_list = np.arange(0.0, total_time, data_collection_interval)
    ax.plot(time_list, env.protein_cumulative_amount)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Protein Amount [bps]')
    ax.set_title('Cumulative Protein Production Plot')
    ax.axis([-5, total_time, -5, env.protein_cumulative_amount[-1]+10])
    ax.grid(True)
    fig.savefig('protein_amount.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    time_list = np.arange(0.0, total_time, dt)
    ax.plot(time_list, env.protein_production)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Protein Amount [bps]')
    ax.set_title('Protein Production Plot')
    ax.axis([0, total_time, -0.5, 4])
    ax.grid(True)
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('protein_production.png')


    fig, ax = plt.subplots(figsize = (10, 5))
    env.plot_velocity_RNAP(ax)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Velocity [bps/s]')
    ax.set_title('RNAP Velocity Plot')
    #ax.axis([-5, 100, -5, 100])
    ax.grid(True)
    fig.savefig('velocity.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    env.plot_velocity_RNAP(ax, start = 0, end = 1)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Velocity [bps/s]')
    ax.set_title('RNAP Velocity Plot')
    ax.axis([-5, 50, -5, 100])
    ax.grid(True)
    fig.savefig('first_velocity.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(range(int(total_time/data_collection_interval)), env.step_processing_time)
    ax.set_xlabel('Stage')
    ax.set_ylabel('Time [s]')
    ax.set_title('Protein Production Plot')
    ax.grid(True)
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('processing_time.png')

    fig, ax = plt.subplots(figsize = (10, 10))
    env.plot_5_and_3_end_amount(ax)
    ax.set_xlabel('Time [min]')
    ax.set_ylabel('Amount')
    ax.set_title('5 End and 3 End Amount Plot')
    ax.grid(True)
    ax.legend()
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('five_and_three.png')

    '''
    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(env.DNA.last_load_location)
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('last_loaded location.png')
    '''

    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(np.arange(0.0, total_time, dt), env.DNA.torq_first_rnap/tau_c)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Torque')
    ax.set_title('Torque of the First RNAP')
    ax.grid(True)
    #ax.set_xlim([0, 250])
    #ax.set_ylim([0,60])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('torq_first_rnap.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(np.arange(0.0, total_time, dt), env.DNA.phi_first_rnap)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('phi')
    ax.set_title('phi of the First RNAP')
    ax.grid(True)
    #ax.set_xlim([0, 250])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('phi_first_rnap.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(np.arange(0.0, total_time, dt), env.DNA.f_n)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('f(n) Value')
    ax.set_title('n_dependence graph')
    ax.grid(True)
    #ax.set_xlim([0, 250])
    #ax.set_ylim([0,250])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('f_n_denpendence.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(np.arange(0.0, total_time, dt), env.DNA.record_promoter_state)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Promoter State')
    ax.set_title('promoter state graph')
    ax.grid(True)
    ax.set_xlim([0, 250])
    ax.set_ylim([0, 2])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('promoter_state.png')

    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(np.arange(0.0, total_time, dt), env.DNA.record_last_r_ref)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Promoter State')
    ax.set_title('promoter state graph')
    ax.grid(True)
    ax.set_xlim([0, 250])
    fig.patch.set_facecolor('xkcd:white')
    fig.savefig('last_r_ref.png')

    end_time = time.perf_counter()
    print(f'Graphing Runtime: {end_time-start_time} s.')
#||---------------------------------------------------------------------||