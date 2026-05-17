try:
    import gwsurrogate as gwsurr
except ModuleNotFoundError as e:
    raise RuntimeError("The gwsurrogate package has failed to load") from e

# VU: temporary until local code merges. 
try:
    import gwsurrogate_eccentric as gwsurr_ecc
except ModuleNotFoundError as e:
    raise RuntimeError("The gwsurrogate_eccentric package has failed to load") from e

import lal
import numpy as np
from gwpy.timeseries import TimeSeries
import astropy.units as u

from lalsimulation.gwsignal.core.waveform import CompactBinaryCoalescenceGenerator
import lalsimulation.gwsignal.core.gw as gw
import lalsimulation as lalsim

# ignore spin magnitude/mass ratio outside training space warnings
import warnings
warnings.filterwarnings("ignore", message=".*Spin")
warnings.filterwarnings("ignore", message=".*Mass ratio")


class NRHybSur3dq8_gwsurr(CompactBinaryCoalescenceGenerator):
    """
    Implements a toy wrapper for NRHybSur3dq8 in the gwsurrogate package
    """

    def __init__(self, **kwargs):
        self.sur = gwsurr.LoadSurrogate("NRHybSur3dq8")
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "aligned-spin",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NRHybSur3dq8",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
            "length": "long"
        }
        return metadata

    def generate_td_modes(self, **parameters):

        # VU: Hacky because gwsignal is being weird about how modes are passed...
        mode_array = parameters.pop("ModeArray")

        # VU: note that this parameter check will convert Mpc to pc
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
        q = m1 / m2  # This is the gwsurrogate convention, q=m1/m2>=1
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
            dist_mpc=dist / 1e6,  # gwsignal converts Mpc to pc
            mode_list = mode_array
        )

        hlm = self._to_gwpy_series(h, times)
        return gw.GravitationalWaveModes(hlm)

    def generate_td_waveform(self, **parameters):
        # VU: added pi/2-phi_ref to match LALSim convention
        theta, phi = (
            parameters["inclination"],
            (np.pi / 2 - parameters["phi_ref"].value) * u.rad,
        )
        hlm = self.generate_td_modes(**parameters)
        hp, hc = hlm(theta, phi)
        hp, hc = TimeSeries(hp, name="hplus"), TimeSeries(hc, name="hcross")
        return hp, hc

    def generate_fd_polarizations_from_td(self, **parameters):
        # VU: inspired by LALSimInspiralGeneratorConditioning.c L486

        # Adjust deltaT depending on sampling rate
        fmax = parameters["f_max"].value
        f_nyquist = fmax
        deltaF = 0
        if "deltaF" in parameters.keys():
            deltaF = parameters["deltaF"].value

        if deltaF != 0:
            n = int(np.round(fmax / deltaF))
            if n & (n - 1):
                chirplen_exp = np.frexp(n)
                f_nyquist = np.ldexp(1, int(chirplen_exp[1])) * deltaF

        deltaT = 0.5 / f_nyquist
        parameters["deltaT"] = deltaT * u.s

        # VU: reduce fmin to make sure tapering doesn't remove signal: [cf. L#1046 in SimInspiral.c]
        m1 = parameters["mass1"].value
        m2 = parameters["mass2"].value
        s1z = parameters["spin1z"].value
        s2z = parameters["spin2z"].value
        fmin = parameters["f22_start"].value
        extra_cycles = 3.0
        extra_time_fraction = 0.1
        m1_kg = m1 * lal.MSUN_SI
        m2_kg = m2 * lal.MSUN_SI

        tchirp = lalsim.SimInspiralChirpTimeBound(fmin, m1_kg, m2_kg, s1z, s2z)
        s = lalsim.SimInspiralFinalBlackHoleSpinBound(s1z, s2z)
        tmerge = lalsim.SimInspiralMergeTimeBound(
            m1_kg, m2_kg
        ) + lalsim.SimInspiralRingdownTimeBound(m1_kg + m2_kg, s)
        textra = extra_cycles / fmin
        fstart = lalsim.SimInspiralChirpStartFrequencyBound(
            (1.0 + extra_time_fraction) * tchirp + tmerge + textra, m1_kg, m2_kg
        )

        # VU: finally update the parameters dict that gets sent to gwsurrogate
        parameters["f22_start"] = fstart * u.Hz
        hp_, hc_ = self.generate_td_waveform(**parameters)

        # VU: reset this in case it gets used somewhere else
        parameters["f22_start"] = fmin * u.Hz

        epoch = lal.LIGOTimeGPS(hp_.times[0].value)

        hp = lal.CreateREAL8TimeSeries(
            "hplus",
            epoch,
            0,
            parameters["deltaT"].value,
            lal.DimensionlessUnit,
            len(hp_),
        )
        hc = lal.CreateREAL8TimeSeries(
            "hcross",
            epoch,
            0,
            parameters["deltaT"].value,
            lal.DimensionlessUnit,
            len(hc_),
        )

        hp.data.data = hp_.value
        hc.data.data = hc_.value

        lalsim.SimInspiralTDConditionStage1(
            hp, hc, extra_time_fraction * tchirp + textra, fmin
        )

        fisco = 1.0 / (
            (6.0**1.5) * lal.PI * (m1_kg + m2_kg) * lal.MTSUN_SI / lal.MSUN_SI
        )
        lalsim.SimInspiralTDConditionStage2(hp, hc, fmin, fisco)

        # Adjust signal duration
        if deltaF == 0:
            chirplen = hp.data.length
            chirplen_exp = np.frexp(chirplen)
            chirplen = int(np.ldexp(1, chirplen_exp[1]))
            deltaF = 1.0 / (chirplen * deltaT)
            parameters["deltaF"] = deltaF

        else:
            chirplen = int(1.0 / (deltaF * deltaT))

        # resize waveforms to the required length
        lal.ResizeREAL8TimeSeries(hp, hp.data.length - chirplen, chirplen)
        lal.ResizeREAL8TimeSeries(hc, hc.data.length - chirplen, chirplen)

        # FFT - Using LAL routines
        hptilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_PLUS",
            hp.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )
        hctilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_CROSS",
            hc.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )

        plan = lal.CreateForwardREAL8FFTPlan(chirplen, 0)
        lal.REAL8TimeFreqFFT(hctilde, hc, plan)
        lal.REAL8TimeFreqFFT(hptilde, hp, plan)

        return hptilde, hctilde

    def _to_gwpy_series(self, modes_dict, times):
        """
        Iterate over the dict and return a dict of gwpy TimeSeries objects
        """
        gwpy_dict = {}
        for ellm, mode in modes_dict.items():
            gwpy_dict[ellm] = TimeSeries(
                mode, times=times, name="h_%i_%i" % (ellm[0], ellm[1])
            )
        return gwpy_dict

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            if isinstance(waveform_dict[key], u.Quantity):
                new_dc[key] = waveform_dict[key].value
            else:
                new_dc[key] = waveform_dict[key]
        return new_dc

class NRHybSur3dq8_CCE_gwsurr(NRHybSur3dq8_gwsurr):
    """
    Implements a toy wrapper for NRHybSur3dq8_CCE in the gwsurrogate package
    """

    def __init__(self, **kwargs):
        self.sur = gwsurr.LoadSurrogate("NRHybSur3dq8_CCE")
        self._update_domains()

class NRSur7dq4_gwsurr(CompactBinaryCoalescenceGenerator):
    """
    Implements a toy wrapper for NRSur7dq4 in the gwsurrogate package
    """

    def __init__(self, **kwargs):

        # super().__init__()
        self.sur = gwsurr.LoadSurrogate("NRSur7dq4")
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "precessing",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NRSur7dq4",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
            "length": "short"
        }
        return metadata

    def generate_td_modes(self, **parameters):
        """
        Generate modes by calling gwsurrogate
        """


        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        lmax = self.waveform_dict.get("ModeArray",None)
        f_start, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]
        f_ref = self.waveform_dict["f22_ref"]

        m1, m2 = self.waveform_dict["mass1"], self.waveform_dict["mass2"]
        chi1 = np.array(
            [
                self.waveform_dict["spin1x"],
                self.waveform_dict["spin1y"],
                self.waveform_dict["spin1z"],
            ]
        )
        chi2 = np.array(
            [
                self.waveform_dict["spin2x"],
                self.waveform_dict["spin2y"],
                self.waveform_dict["spin2z"],
            ]
        )
        dist = self.waveform_dict["distance"]
        q = m1 / m2  # This is the gwsurrogate convention, q=m1/m2>=1
        if q < 1.0:
            raise Exception("m2 should not be bigger than m1!")

        times, h, dyn = self.sur(
            q,
            chi1,
            chi2,
            dt=dt,
            f_low=f_start,
            f_ref=f_ref,
            units="mks",  # Output in SI units
            M=m1 + m2,  # In solar masses
            dist_mpc=dist / 1e6,  # In Mpc
            ellMax = lmax
        )

        hlm = self._to_gwpy_series(h, times)
        return gw.GravitationalWaveModes(hlm)

    def generate_td_waveform(self, **parameters):
        """
        Generate plus and cross polarizations from the modes
        """
        theta, phi = (
            parameters["inclination"],
            (np.pi / 2 - parameters["phi_ref"].value) * u.rad,
        )
        hlm = self.generate_td_modes(**parameters)
        hp, hc = hlm(theta, phi)
        hp, hc = TimeSeries(hp, name="hplus"), TimeSeries(hc, name="hcross")
        return hp, hc

    def generate_fd_polarizations_from_td(self, **parameters):
        """
        Generate frequency domain plus and cross polarizations for parameters
        VU: Inspired by LALSimInspiralGeneratorConditioning.c L486
        since this is a short waveform
        """

        # Adjust deltaT depending on sampling rate
        fmax = parameters["f_max"].value
        f_nyquist = fmax
        deltaF = 0
        if "deltaF" in parameters.keys():
            deltaF = parameters["deltaF"].value

        if deltaF != 0:
            n = int(np.round(fmax / deltaF))
            if n & (n - 1):
                chirplen_exp = np.frexp(n)
                f_nyquist = np.ldexp(1, int(chirplen_exp[1])) * deltaF

        deltaT = 0.5 / f_nyquist
        parameters["deltaT"] = deltaT * u.s

        hp_, hc_ = self.generate_td_waveform(**parameters)

        epoch = lal.LIGOTimeGPS(hp_.times[0].value)

        hp = lal.CreateREAL8TimeSeries(
            "hplus",
            epoch,
            0,
            parameters["deltaT"].value,
            lal.DimensionlessUnit,
            len(hp_),
        )
        hc = lal.CreateREAL8TimeSeries(
            "hcross",
            epoch,
            0,
            parameters["deltaT"].value,
            lal.DimensionlessUnit,
            len(hc_),
        )

        hp.data.data = hp_.value
        hc.data.data = hc_.value

        # conditioning/tapering is done differently since this is a short waveform
        # [cf. L#44 in LALSimInspiralGeneratorConditioning.c]
        taper = True
        lalsim.SimInspiralREAL8WaveTaper(hp.data, taper)
        lalsim.SimInspiralREAL8WaveTaper(hc.data, taper)

        # Adjust signal duration
        if deltaF == 0:
            chirplen = hp.data.length
            chirplen_exp = np.frexp(chirplen)
            chirplen = int(np.ldexp(1, chirplen_exp[1]))
            deltaF = 1.0 / (chirplen * deltaT)
            parameters["deltaF"] = deltaF

        else:
            chirplen = int(1.0 / (deltaF * deltaT))

        # resize waveforms to the required length
        lal.ResizeREAL8TimeSeries(hp, hp.data.length - chirplen, chirplen)
        lal.ResizeREAL8TimeSeries(hc, hc.data.length - chirplen, chirplen)

        # FFT - Using LAL routines
        hptilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_PLUS",
            hp.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )
        hctilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_CROSS",
            hc.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )

        plan = lal.CreateForwardREAL8FFTPlan(chirplen, 0)
        lal.REAL8TimeFreqFFT(hctilde, hc, plan)
        lal.REAL8TimeFreqFFT(hptilde, hp, plan)

        return hptilde, hctilde

    def _to_gwpy_series(self, modes_dict, times):
        """
        Iterate over the dict and return a dict of gwpy TimeSeries objects
        """
        gwpy_dict = {}
        for ellm, mode in modes_dict.items():
            gwpy_dict[ellm] = TimeSeries(
                mode, times=times, name="h_%i_%i" % (ellm[0], ellm[1])
            )
        return gwpy_dict

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            if isinstance(waveform_dict[key], u.Quantity):
                new_dc[key] = waveform_dict[key].value
            else:
                new_dc[key] = waveform_dict[key]
        return new_dc

class NR_hdf5_gwsurr(CompactBinaryCoalescenceGenerator):
    """
    Implements a wrapper for the NR_hdf5 model called from LALSim
    """

    def __init__(self, **kwargs):

        super().__init__()
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "precessing",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NR_hdf5",
            "implementation": "C++",
            "conditioning_routines": "gwsignal",
            "length": "short"
        }
        return metadata

    def generate_td_modes(self, **parameters):
        """
        Generate modes by calling LALSim for NR td waveform
        """
        return 

    def generate_td_waveform(self, **parameters):
        """
        Generate plus and cross polarizations from SimInspiralTD
        """

        mode_array = parameters.pop("ModeArray")
        nr_data_file = parameters.pop('NumRelData')

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        f_start, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]
        f_ref = self.waveform_dict["f22_ref"]

        mass1 = self.waveform_dict["mass1"]
        mass2 = self.waveform_dict["mass2"]

        spin1x = self.waveform_dict["spin1x"]
        spin1y = self.waveform_dict["spin1y"]
        spin1z = self.waveform_dict["spin1z"]
        spin2x = self.waveform_dict["spin2x"]
        spin2y = self.waveform_dict["spin2y"]
        spin2z = self.waveform_dict["spin2z"]

        luminosity_distance= self.waveform_dict["distance"] # stripped from u.pc

        theta, phi = (
            self.waveform_dict["inclination"],
            np.pi / 2 - self.waveform_dict["phi_ref"],
        )

        approximant = lalsim.NR_hdf5

        # Setup LAL parameters
        LALparams = lal.CreateDict()
        lalsim.SimInspiralWaveformParamsInsertNumRelData(LALparams, nr_data_file)

        if mode_array is not None:
            mode_array_lal = lalsim.SimInspiralCreateModeArray()
            for mode in mode_array:
                lalsim.SimInspiralModeArrayActivateMode(mode_array_lal, mode[0], mode[1])
            lalsim.SimInspiralWaveformParamsInsertModeArray(LALparams, mode_array_lal)

        # Get waveform polarizations from LALSim
        hp_lal, hc_lal = lalsim.SimInspiralChooseTDWaveform(
            mass1 * lal.MSUN_SI,
            mass2 * lal.MSUN_SI,
            spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z,
            luminosity_distance * lal.PC_SI,
            theta,
            phi,
            0.,
            0.,
            0.,
            dt,
            f_start,
            f_ref,
            LALparams,
            approximant
        )

        hp_data = hp_lal.data.data
        hc_data = hc_lal.data.data

        hp, hc = TimeSeries(hp_data, name="hplus"), TimeSeries(hc_data, name="hcross")
        return hp, hc

    def generate_fd_polarizations_from_td(self, **parameters):
        """
        """

        # Adjust deltaT depending on sampling rate
        fmax = parameters["f_max"].value
        f_nyquist = fmax
        deltaF = 0
        if "deltaF" in parameters.keys():
            deltaF = parameters["deltaF"].value

        if deltaF != 0:
            n = int(np.round(fmax / deltaF))
            if n & (n - 1):
                chirplen_exp = np.frexp(n)
                f_nyquist = np.ldexp(1, int(chirplen_exp[1])) * deltaF

        deltaT = 0.5 / f_nyquist
        parameters["deltaT"] = deltaT * u.s

        hp_, hc_ = self.generate_td_waveform(**parameters)

        epoch = lal.LIGOTimeGPS(hp_.times[0].value)

        hp = lal.CreateREAL8TimeSeries(
            "hplus",
            epoch,
            0,
            parameters["deltaT"].value,
            lal.DimensionlessUnit,
            len(hp_),
        )
        hc = lal.CreateREAL8TimeSeries(
            "hcross",
            epoch,
            0,
            parameters["deltaT"].value,
            lal.DimensionlessUnit,
            len(hc_),
        )

        hp.data.data = hp_.value
        hc.data.data = hc_.value

        # conditioning/tapering is done differently since this is a short waveform
        # [cf. L#44 in LALSimInspiralGeneratorConditioning.c]
        taper = True
        lalsim.SimInspiralREAL8WaveTaper(hp.data, taper)
        lalsim.SimInspiralREAL8WaveTaper(hc.data, taper)

        # Adjust signal duration
        if deltaF == 0:
            chirplen = hp.data.length
            chirplen_exp = np.frexp(chirplen)
            chirplen = int(np.ldexp(1, chirplen_exp[1]))
            deltaF = 1.0 / (chirplen * deltaT)
            parameters["deltaF"] = deltaF

        else:
            chirplen = int(1.0 / (deltaF * deltaT))

        # resize waveforms to the required length
        lal.ResizeREAL8TimeSeries(hp, hp.data.length - chirplen, chirplen)
        lal.ResizeREAL8TimeSeries(hc, hc.data.length - chirplen, chirplen)

        # FFT - Using LAL routines
        hptilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_PLUS",
            hp.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )
        hctilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_CROSS",
            hc.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )

        plan = lal.CreateForwardREAL8FFTPlan(chirplen, 0)
        lal.REAL8TimeFreqFFT(hctilde, hc, plan)
        lal.REAL8TimeFreqFFT(hptilde, hp, plan)

        return hptilde, hctilde

    def _to_gwpy_series(self, modes_dict, times):
        """
        Iterate over the dict and return a dict of gwpy TimeSeries objects
        """
        gwpy_dict = {}
        for ellm, mode in modes_dict.items():
            gwpy_dict[ellm] = TimeSeries(
                mode, times=times, name="h_%i_%i" % (ellm[0], ellm[1])
            )
        return gwpy_dict

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            if isinstance(waveform_dict[key], u.Quantity):
                new_dc[key] = waveform_dict[key].value
            else:
                new_dc[key] = waveform_dict[key]
        return new_dc

class NRSurE_NoSpin_gwsurr(NRSur7dq4_gwsurr):
    """
    Implements a toy wrapper for NRSur7dq4 in the gwsurrogate package
    """

    def __init__(self, **kwargs):

        # super().__init__()

        self.sur = gwsurr_ecc.LoadSurrogate('NRSurE_q4NoSpin_22')
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "precessing",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NRSurE_??",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
            "length": "short"
        }
        return metadata

    def generate_td_modes(self, **parameters):
        """
        Generate modes by calling gwsurrogate
        """


        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        f_start, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]

        if f_start!=0:
            raise Exception('KILL NRSur_E does not currently support nonzero starting frequencies')

        # VU: f_ref isn't actually passed, this might create problems during PE
        f_ref = self.waveform_dict["f22_ref"]

        m1, m2 = self.waveform_dict["mass1"], self.waveform_dict["mass2"]
        chi1 = np.array([0.,0.,0.])
        chi2 = np.array([0.,0.,0.])
        e = self.waveform_dict['eccentricity']
        meanano = self.waveform_dict['meanPerAno']
        dist = self.waveform_dict["distance"]
        q = m1 / m2  # This is the gwsurrogate convention, q=m1/m2>=1
        if q < 1.0:
            raise Exception("m2 should not be bigger than m1!")

        times, h, dyn = self.sur(
            q,
            chi1,
            chi2,
            e=e,
            meanano=meanano,
            dt=dt,
            f_low=f_start,
            units="mks",  # Output in SI units
            M=m1 + m2,  # In solar masses
            dist_mpc=dist / 1e6,  # In Mpc
        )

        hlm = self._to_gwpy_series(h, times)
        return gw.GravitationalWaveModes(hlm)

class NRSur7dq4v2_gwsurr(NRSur7dq4_gwsurr):
    """
    Implements a toy wrapper for NRSur7dq4v2 in the gwsurrogate package
    """

    def __init__(self, **kwargs):

        # super().__init__()
        self.sur = gwsurr.LoadSurrogate("NRSur7dq4v2")
        self._update_domains()
    
    @property
    def metadata(self):
        # Override the metadata to ensure correct pipeline identification
        metadata = super().metadata
        metadata["approximant"] = "NRSur7dq4v2"
        return metadata
