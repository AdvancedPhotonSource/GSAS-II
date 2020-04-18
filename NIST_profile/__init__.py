from __future__ import absolute_import
from .profile_functions_class import FP_profile, best_irfft, best_rfft
from .atan_windowed_FP_profile import FP_atan_windowed_convolver

class FP_windowed(FP_atan_windowed_convolver, FP_profile):
    pass


