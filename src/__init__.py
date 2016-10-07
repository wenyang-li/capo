
#import gb_weather, pfb, pspec, dspec, red, fringe, omni, frf_conv, uCal, miriad
#import arp, jcp, dfm, dcj, zsa, ctc, wyl

import pfb, pspec, dspec, red, fringe, miriad
import fringe as frf_conv # for backward compatibility
import oqe
import warnings
try: import omni, uCal
except(ImportError,NameError):
    warnings.warn("Warning: omnical not installed, not importing capo.omni")
import arp, jcp, dfm, dcj, zsa, ctc, wyl

