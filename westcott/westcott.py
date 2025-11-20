from .tables import *
from .user import *
from .westcott_gfactors import *

class Westcott(gFactors):
    __doc__="""Class to handle Westcott g-factor calculations."""

    def __init__(self):
        CrossSectionData.__init__(self)
        NeutronFlux.__init__(self)
        ResonanceData.__init__(self)
        UserSpectrum.__init__(self)
        Kinematics.__init__(self)
        Irregularity.__init__(self)
        gFactors.__init__(self)
