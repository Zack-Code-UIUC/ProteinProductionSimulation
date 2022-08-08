"""
=======
mrna.py
=======

This file represent the messanger RNA, abbreviated as mRNA. The mRNA is produced from the transcription of DNA by the
RNAP. Ribosomes will attach to mRNA to conduct translation to synthesize proteins.


"""

from ..interface import Entity


class mRNA(Entity):
    """
    Represent one mRNA. The mRNA class also manage all the ribosome and protein related action.
    """