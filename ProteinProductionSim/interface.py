"""
============
interface.py
============

This file lists and define all the interfaces that is used in this package. There are four types of interfaces, being
    * Controller
    * Data Container
    * Environment
    * Entity

---------------------------------
How each interface is structured?
---------------------------------


init()
^^^^^^^

Well, all class must have __init__ method, for sure. In addition, each interface also must have init method.

The init method is necessary, because the user might modify the Controller instance after the initiation of the
instance. Therefore, each class must initiate its relevant quantities in its init method to counter any random change
to its quantities by the user.

step()
^^^^^^

For the Controller, Environment and Entity class, they all have the step method. The step method essentially represent a
single stepping. In principle, the Controller.step() will call the Environment.step(). And, the Environment.step() will
call the Entity.step().

parent
^^^^^^

For the DataContainer, Environment and Entity class, they all have the class variable "parent". This variable points
back to its parent instance, in order to use callback method.

_call_back()
^^^^^^^^^^^^

Except tbe Controller class, the _call_back() method is intended to be called by its children instances for information
transfer.

-------------------------------------
How each class relates to each other?
-------------------------------------

            Controller
                ||
                || Record information using DataContainer
                ||
            Environment
                ||
               /||\
            // // \\ \\
            Many Entities


The Controller class contains all other classes. It uses the DataContainer class to log and store information. The
Controller class also controls the Environment class, which contains all the entities instance. This is a tree.

There are reasons behind such tree structure, coming rom perspectives of abstraction, extendability, and malleability.
Originally, the classes Controller, DataContainer and Environment are all concentrated into one class, but that prove to
too much duty in one hand. The code for this old class grow very big and annoying to look at. :)

Now, the Controller acts as the switch to the start the simulation. The DataContainer is used for data. The Environment
represents the "Environment" or the simulated contents. The entity represents the all the simulated entities, like a
DNA strand or a mRNA.
"""


class Controller:
    """
    The Controller interface is intended to contain all the DataContainer instance and used to start the simulation.
    """
    def __init__(self):
        pass

    def init(self):
        pass

    def start(self):
        pass

    def call_back(self, option, data):
        pass


class DataContainer:
    """
    The DataContainer interface is used to log and store data.
    """

    def __init__(self, parent, **kwargs):
        self.parent = parent
        pass

    def init(self, **kwargs):
        pass

    def log(self, **kwargs):
        pass

    def get(self, **kwargs):
        pass

    def plot(self, **kwargs):
        pass


class Environment:
    """
    The Environment interface is used to represent the simulated environment and contains all the entities instance.
    """

    def __init__(self, parent):
        self.parent = parent
        pass

    def step(self, **kwargs):
        pass

    def init(self):
        pass

    def call_back(self, option, data):
        pass


class Entity:
    """
    The Entity interface is used to represent the simulated entities.
    """

    def __init__(self, parent, **kwargs):
        self.parent = parent
        pass

    def init(self):
        pass

    def step(self, **kwargs):
        pass

    def call_back(self, option, data=None):
        pass
