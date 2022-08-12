

def make_directory_now(direc, name):
    '''
    Create a directory based on the target directory and inputed name. A date and time string will be attached to the back of the name.

    Parameters:
    direc(string): target directory
    name(string) : directory name that need to be created.

    Return:
    string: the path of the create directory
    '''
    # get current working directory
    cwd = os.getcwd()
    # construct the date and time string
    now = datetime.now()
    date_string = now.strftime("%Y_%m_%d_%H_%M")
    # make the directory name
    directory_name = nam e +'_ ' +date_string
    # construct the whole path
    path_directory = cw d +'/ ' +dire c +'/ ' +directory_name
    try:
        makedirs(path_directory)
    except FileExistsError:
        pass
    return path_directory