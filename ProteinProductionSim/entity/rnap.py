"""
================================
rnap.py
================================

This file represents the RNA Polymerase, abbreviated as RNAP. The RNAP is attached to
the DNA strand to conduct the transcription of messanger RNA.


"""
from ..interface import Entity
from ..helper.random_generator import binary_generator, exponential_generator, stepwise_exponential_generator
from ..helper.loading_list import LoadingList
from datacontainer.ribo_container import RIBOContainer
from ..variables import length, kRiboLoading, ribo_loading_interval, m1, m2, t_crit, \
    initiation_nt, scaling, dt


class RNAP(Entity):
    """Represent the RNA polymerase.

    Parameters
    ---------
    dna : DNAStrand
        the reference of the parent DNA strand
    initial_t : float
        the initial time when this RNAP is attached  or loaded to the DNA strand.
    pause_profile : str, optional
        the site-specific pausing pattern that is used (default is 'flat')
    ribo_loading_profile : str, optional
        this loading pattern for the ribosomes (default is "stochastic")
    degradation_profile : str, optional
        the degradation time or loading interval pattern that is used (default is "exponential")


    Attributes
    ----------
    parent : DNAStrand
        the reference to the parent DNAStrand
    initial_t : int
        the initial time when this RNAP is attached  or loaded to the DNA strand
    position : float
        the position of the RNAP on the DNA
    attached : bool
        this shows if the RNAP is attached to the DNA
    passed_site_1 : bool
        this shows if the RNAP has passed the pausing site 1
    passed_site_2 : bool
    passing_1 : bool
    passing_2 : bool
    initiated : bool
    degradation_profile : str
    degrading : bool
    degraded : bool
    t_degrade : int
    loading_list : LoadingList
    RIBO_LIST : List




    """
    pauseProb = 0.8

    def __init__(self, dna: Entity, initial_t, pause_profile: str = "flat", ribo_loading_profile: str = "stochastic",
                 degradation_profile: str = "exponential", protein_production_off: bool = False):
        super().__init__(dna)
        self.parent = dna  # this store the reference to its mother DNA, so that callback method can be used.
        self.initial_t = initial_t  # means the initial time when the RNAP attached.
        self.position = 0  # means position.
        self.attached = True  # indicated if the RNAP is attached to the DNA. Will detached if reach the end.

        # Site-Pausing
        self.passed_site_1 = False  # passed_site_1 is indicating whether the RNAP has passed the pausing site.
        self.passed_site_2 = False  # similar to above. This two variable is also used to bypass mechanisms.
        self.passing_1 = False  # if the RNAP is passing pausing site 1.
        self.passing_2 = False  # if the RNAP is passing pausing site 2.
        match pause_profile:
            case "flat":
                self.passed_site_1 = True
                self.passed_site_2 = True
            case "OnepauseAbs":
                self.passed_site_1 = binary_generator(1-RNAP.pauseProb)
                self.passed_site_2 = True
            case "TwopauseAbs":
                self.passed_site_1 = binary_generator(1-RNAP.pauseProb)
                self.passed_site_2 = binary_generator(1-RNAP.pauseProb)

        # mRNA degradation
        self.initiated = False  # initiated indicates if the length has passed the size required for initiation (33nts)

        # mRNA degradation
        self.degradation_profile = degradation_profile
        self.degrading = False
        self.degraded = False
        match self.degradation_profile:
            case "exponential":
                self.t_degrade = scaling(exponential_generator(1 / ribo_loading_interval))
            case "stepwise exponential":
                self.t_degrade = scaling(stepwise_exponential_generator(m1, m2, t_crit))

        # Loading of Ribosomes
        match ribo_loading_profile:
            case "uniform":
                self.loading_list = LoadingList(self, scaling(self.t_degrade), kRiboLoading*dt, if_stochastic=False,
                                                if_bursty=False)
            case "stochastic":
                self.loading_list = LoadingList(self, scaling(self.t_degrade), kRiboLoading*dt, if_stochastic=True,
                                                if_bursty=False)

        # sometimes we do not want to activate the protein production, then we just dump the whole loading list.
        if protein_production_off:
            self.loading_list.dump()

        # we use the DataContainer RIBOContainer to both store and manage the Ribosomes
        # the class RIBOContainer will contain various class functions to helpe us with RIBO-related business
        self.RIBO_LIST = RIBOContainer(self)

    def init(self):
        """
        This init() method is empty.
        """
        pass

    def step(self, time_index: int, pace: float) -> int:
        """
        This method step the RNAP and its children mRNA instance forward.

        Parameters
        ----------
        time_index : int
            this represents the current time point
        pace : float
            this represents how much distance has the RNAP moved on the DNA.

        RETURNS
        -------
        int : the amount of proteins produced
        """
        # REASON: check whether the mRNA is degraded, if so, then skip this instance completely.
        if self.degraded:
            return 0

        # REASON: increment the RNAP position by the amount pace.
        self.position += pace

        # REASON: check if the RNAP is detached, if so, then set the RNAP position to a ridiculously far place.
        #         use callback method of the parent DNA to increment the detached amount.
        if self.attached and (self.position >= length):
            self.attached = False
            self.position = length + 1000
            self.parent.call_back(option="attached", data=-1)
            self.parent.call_back(option="detached")

        # REASON: check for degradation initiation. check for initiation, which allow the loading of Ribosome.
        #         check for loading of Ribosome on the mRNA.
        if not self.degrading:
            # REASON: we immediately check for degradation.
            #         if the time has exceeded the degradation time, then the degrading state will be mark True.
            if time_index >= self.t_degrade + self.initial_t:
                self.degrading = True
                self.parent.call_back(option="degrading")

            # REASON: check for self.initiated,
            #         if the RNAP instance is not initiated, and it satisfies to requirement to be initiated
            #         then we initiate it :)
            if not self.initiated and self.position >= initiation_nt:
                if self.position >= initiation_nt:
                    self.initiated = True

            # REASON: if the mRNA is initiated and the loading list is not empty, then check for potential Loading.
            if self.initiated and self.loading_list.if_empty():
                # REASON: check with the loading_list if we can load a ribosome right not.
                if self.loading_list.if_can_load(time_index, self.initial_t):
                    # REASON: check for potential hindrance
                    if self.RIBO_LIST.if_clear_at_start():
                        self.RIBO_LIST.load_one()

        # REASON: we let the RIBO_LIST handle the stepping of the Ribosome
        #         we can trust it to check for hindrance and various matters
        prot = self.RIBO_LIST.step(time_index, self.position)

        # REASON: check for complete degradation. if completely degraded, then set self.degraded to True
        #         first the RNAP has to be detached.
        #         second the mRNA has to be degrading
        #         third all the ribosome has to be already detached from the mRNA
        if not self.attached and self.degrading and self.RIBO_LIST.if_all_detached():
            self.degraded = True
            self.parent.call_back(option="degraded")
        # return protein production
        return prot

    def call_back(self, option, data):
        """
        This call_back() method is empty and does not perform any task.
        """
        pass
