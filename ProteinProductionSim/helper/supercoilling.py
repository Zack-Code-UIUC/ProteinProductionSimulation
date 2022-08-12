from variables import gamma


def s(r):
    return gamma * r


def n_dependence_cubic_3(x):
    # return 1+0.778753*(x-1)+3.3249*(x-1)**2+0.379478*(x-3)**3
    return 1+0.778753*(x-1)+3.3249*(x-1)**2+0.379478*(x-1)**3
