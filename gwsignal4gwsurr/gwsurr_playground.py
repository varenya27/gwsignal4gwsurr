"""
Testing grounds for models that are either being developed or 
just made for playing around/debugging purposes
"""

import gwsurrogate as gwsurr
from .gwsurr import NRHybSur3dq8_gwsurr,NRSur7dq4_gwsurr
import numpy as np
import lal
import lalsimulation as lalsim
import lalsimulation.gwsignal.core.gw as gw

# ignore spin magnitude/mass ratio outside training space warnings
import warnings
warnings.filterwarnings("ignore", message=".*Spin")
warnings.filterwarnings("ignore", message=".*Mass ratio")

class NRSur3dq8_Lev2_varenya_gwsurr(NRHybSur3dq8_gwsurr):
    def __init__(self, **kwargs):
        self.sur = gwsurr.LoadSurrogate("NRSur3dq8_Lev2_varenya")
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "aligned-spin",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NRSurr",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
        }
        return metadata

    def generate_td_modes(self, **parameters):

        extra_args = parameters.pop('extra_args')
        noisy = extra_args.pop("noisy")
        noise_level = extra_args.pop("noise_level")

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        fstart, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]
        f_ref = self.waveform_dict["f22_ref"]

        m1, m2 = self.waveform_dict["mass1"], self.waveform_dict["mass2"]
        s1z = self.waveform_dict["spin1z"]
        s2z = self.waveform_dict["spin2z"]
        chi1 = np.array(
            [
                0.0,
                0.0,
                s1z,
            ]
        )
        chi2 = np.array(
            [
                0.0,
                0.0,
                s2z,
            ]
        )
        dist = self.waveform_dict["distance"]
        q = m1 / m2
        if q < 1.0:
            raise Exception("m2 should not be bigger than m1!")

        times, h, dyn = self.sur(
            q,
            chi1,
            chi2,
            dt=dt,
            f_low=fstart,
            f_ref=f_ref,
            units="mks",  # Output in SI units
            M=m1 + m2,  # In solar masses
            dist_mpc=dist / 1e6,  # In Mpc
            mode_list = [(2,2)]
        )

        #=============== add some phase noise ===============#
        if noisy:
            delta_phi = noise_level * np.random.randn(len(h[(2, 2)]))
            for ellm, h_array in h.items():
                ell, m = ellm
                h[ellm] = h_array * np.exp(1j * m * delta_phi)
        #====================================================#

        hlm = self._to_gwpy_series(h, times)
        return gw.GravitationalWaveModes(hlm)


class NRSur3dq8_Lev3_varenya_gwsurr(NRSur3dq8_Lev2_varenya_gwsurr):
    def __init__(self, **kwargs):
        self.sur = gwsurr.LoadSurrogate("NRSur3dq8_Lev3_varenya")
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "aligned-spin",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NRSurr",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
        }
        return metadata

class NRSur7dq4_LALSim_gwsurr(NRSur7dq4_gwsurr):
    def __init__(self, **kwargs):
        self.sur = gwsurr.LoadSurrogate("NRSur7dq4")
        self._update_domains()
        super().__init__()

    def generate_fd_polarizations(self, **parameters):
        parameters = self._strip_units(parameters)
        f_start = parameters["f22_start"]
        f_ref = parameters["f22_ref"]

        mass1, mass2 = parameters["mass1"], parameters["mass2"]
        spin1x = parameters["spin1x"]
        spin1y = parameters["spin1y"]
        spin1z = parameters["spin1z"]
        spin2x = parameters["spin2x"]
        spin2y = parameters["spin2y"]
        spin2z = parameters["spin2z"]
        distance = parameters["distance"]

        hp_gwsignal, hc_gwsignal = lalsim.SimInspiralFD(
            mass1 * lal.MSUN_SI,
            mass2 * lal.MSUN_SI,
            spin1x,
            spin1y,
            spin1z,
            spin2x,
            spin2y,
            spin2z,
            distance * 1.0e6 * lal.PC_SI,
            parameters['inclination'],
            parameters['phi_ref'],
            0.0,
            0.0,
            0.0,
            parameters['deltaF'],
            f_start,
            parameters['f_max'],
            f_ref,  # f_low, f_max, f_ref
            lal.CreateDict(),
            lalsim.GetApproximantFromString("NRSur7dq4")
        )
        return hp_gwsignal, hc_gwsignal 

