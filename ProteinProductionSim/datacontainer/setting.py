"""
==========
setting.py
==========

This class defines the setting class, which store the setting.

"""

from ..interface import DataContainer


class Setting(DataContainer):
    """
    This stores the setting of a run. This is used to passed all the setting around.
    """
    def __init__(self,
                 total_time=300.0,
                 dt = 1/30,
                 data_collection_interval=1/10,
                 data_smoothing_interval=5,
                 length=3072,
                 t_on=7.8,
                 stop=90,
                 pause_profile="flat",
                 loading_profile="stochastic",
                 degradation_profile="exponential"
                 ):
        self.parent = None
        pass

    def _init(self):
        pass

    def log(self):
        pass

    def get(self):
        pass

    def plot(self):
        pass
