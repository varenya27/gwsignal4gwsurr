"""
Testing grounds for models that are either being developed or 
just made for playing around/debugging purposes
"""

import gwsurrogate as gwsurr
import gwtools
from .gwsurr import NRHybSur3dq8_gwsurr,NRSur7dq4_gwsurr
import numpy as np
import lal
import lalsimulation as lalsim
import lalsimulation.gwsignal.core.gw as gw
from lalsimulation.gwsignal.core.waveform import CompactBinaryCoalescenceGenerator
import astropy.units as u
import warnings
from gwpy.timeseries import TimeSeries

try: 
    import sxs
    from sxs.waveforms.format_handlers.lvc import to_lvc_conventions
except ModuleNotFoundError:
    warnings.warn("Could not find the sxs package, SXS-NR injections will not work")

# ignore spin magnitude/mass ratio outside training space warnings
warnings.filterwarnings("ignore", message=".*Spin")
warnings.filterwarnings("ignore", message=".*Mass ratio")
# warnings.filterwarnings("ignore", message=".*Params")

class NRSur3dq8_Lev2_varenya_gwsurr(NRSur7dq4_gwsurr):
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

        extra_args = parameters.pop('extra_args',None)
        noisy = extra_args.pop("noisy")
        noise_level = extra_args.pop("noise_level")
        noise_model = extra_args.pop("noise_model")

        mode_array = parameters.pop('ModeArray')

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
            mode_list = mode_array
        )

        #=============== add some phase noise ===============#
        if noisy:
            if noise_model == 'gaussian':
                delta_phi = noise_level * np.random.randn(len(h[(2, 2)]))
                for ellm, h_array in h.items():
                    m = ellm[1]
                    h[ellm] = h_array * np.exp(1j * m * delta_phi)
            elif noise_model == 'error_surrogate_phase':
                raise Exception('Not yet implemented')
                    
            else:
                raise Exception('Unrecognized error model')
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

class NRHybSur3dq8_short_gwsurr(NRSur7dq4_gwsurr):
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
            "approximant": "NRSurr",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
        }
        return metadata

    def generate_td_modes(self, **parameters):

        mode_array = parameters.pop('ModeArray')

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
            mode_list = mode_array
        )

        hlm = self._to_gwpy_series(h, times)
        return gw.GravitationalWaveModes(hlm)

# !!!! CODE WORK AHEAD !!!!
class NR_SXS(NRSur7dq4_gwsurr):
    """
    Implements a toy wrapper for NR waveforms using the sxs package
    """

    def __init__(self, **kwargs):
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "precessing",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NR_SXS",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
            "length": "short"
        }
        return metadata

    def generate_td_modes(self, **parameters):
        """
        Generate modes by calling sxs
        """

        SXS_ID = parameters.pop('SXS_ID')
        ModeArray = parameters.pop('ModeArray',None)
        remove_memory = parameters.pop('remove_memory', True)

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)

        fstart, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]
        print('DBUG gwsign sxs got fstart=',fstart)
        f_ref = self.waveform_dict["f22_ref"]

        lmax = self.waveform_dict.get("lmax",None)

        dist_mpc = self.waveform_dict['distance']/1e6
        Mtot = self.waveform_dict['mass1']+self.waveform_dict['mass2']

        # scale factors for conversion to/from SI units 
        amp_scale = Mtot*gwtools.Msuninsec*gwtools.c/(1e6*dist_mpc*gwtools.PC_SI)
        t_scale = gwtools.Msuninsec * Mtot

        # convert GW frequency in Hz to orbital frequency in cycles/M
        fstart_M = fstart*t_scale/2
        f_ref_M = f_ref*t_scale/2
        dt_M = dt /t_scale

        # sxs function call to get dimless strain
        try:
            sxs_bbh = sxs.load('SXS:BBH:%04d'%SXS_ID)

            print('INFO SXS waveform diagnostics:')
            print(f'    SXS ID : SXS:BBH:{SXS_ID:04d} ')
            print(f'    Lev    : {getattr(sxs_bbh, "Lev", "N/A")} ')
            print(f'    ellmax : {lmax} ({getattr(getattr(sxs_bbh, "h", None), "ell_max", "N/A")} available)')
            print('INFO SXS reference parameters:')
            print(f'    mass1  : {sxs_bbh.metadata.get("reference_mass1", "N/A")}')
            print(f'    mass2  : {sxs_bbh.metadata.get("reference_mass2", "N/A")}')
            print(f'    q      : {sxs_bbh.metadata.get("reference_mass_ratio", "N/A")}')
            print(f'    chi1   : {sxs_bbh.metadata.get("reference_dimensionless_spin1", "N/A")}')
            print(f'    chi2   : {sxs_bbh.metadata.get("reference_dimensionless_spin2", "N/A")}')
            print(f'    ecc    : {sxs_bbh.metadata.get("reference_eccentricity", "N/A")}')
            print(f'    ell    : {sxs_bbh.metadata.get("reference_mean_anomaly", "N/A")}')
        except Exception:
            raise Exception('KILL Could not load simulation! Try checking the ID?')

        if not remove_memory:
            raise NotImplementedError

        else:
            print('INFO waveform will not have memory')

            strain = sxs_bbh.load_waveform(
                *sxs_bbh.strain_path,
                transform_to_inertial=True,
            )

            strain_no_memory = strain.remove_memory(sxs_bbh.metadata.relaxation_time)
            times, h, dyn = to_lvc_conventions(
                strain_no_memory.to_corotating_frame(),
                sxs_bbh.horizons,
                f_ref=f_ref_M,
                dt=dt_M,
                f_low=fstart_M,
                ell_max=lmax,
                phi_ref=None,
                inclination=None,
            )

        # sxs does not trim the waveform at f_low/t_low
        # so we take matters into our own hands
        start_ind = np.argmin(np.abs(times-dyn['t_low']))

        # convert to unitful qty using Mtot and d_lum
        times *= t_scale
        h_dict = {}
        for k,v in h.items():
            l,m = k
            if l<2: 
                continue
            if m%2==0:
                h_dict[k] = -v[start_ind:]*amp_scale
            else:
                h_dict[k] = v[start_ind:]*amp_scale

        hlm = self._to_gwpy_series(h_dict, times[start_ind:])
        return gw.GravitationalWaveModes(hlm)

    def generate_td_waveform(self, **parameters):
        theta, phi = (
            parameters["inclination"],
            # parameters["phi_ref"]
            (np.pi / 2 - parameters["phi_ref"].value) * u.rad,
        )
        hlm = self.generate_td_modes(**parameters)
        hp, hc = hlm(theta, phi)
        hp, hc = TimeSeries(hp, name="hplus"), TimeSeries(hc, name="hcross")
        return hp, hc

class NR_SimulationAnnex(NRSur7dq4_gwsurr):
    """
    Implements a toy wrapper for NR waveforms using the sxs package
    """

    def __init__(self, **kwargs):
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "precessing",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "NR_SXS",
            "implementation": "Python",
            "conditioning_routines": "gwsignal",
            "length": "short"
        }
        return metadata

    def generate_td_modes(self, **parameters):
        """
        Generate modes by calling sxs
        """

        SXS_ID = parameters.pop('SXS_ID')
        ModeArray = parameters.pop('ModeArray',None)
        remove_memory = parameters.pop('remove_memory', True)

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)

        fstart, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]
        print('DBUG gwsign sxs got fstart=',fstart)
        f_ref = self.waveform_dict["f22_ref"]

        lmax = self.waveform_dict.get("lmax",None)

        dist_mpc = self.waveform_dict['distance']/1e6
        Mtot = self.waveform_dict['mass1']+self.waveform_dict['mass2']

        # scale factors for conversion to/from SI units 
        amp_scale = Mtot*gwtools.Msuninsec*gwtools.c/(1e6*dist_mpc*gwtools.PC_SI)
        t_scale = gwtools.Msuninsec * Mtot

        # convert GW frequency in Hz to orbital frequency in cycles/M
        fstart_M = fstart*t_scale/2
        f_ref_M = f_ref*t_scale/2
        dt_M = dt /t_scale

        return gw.GravitationalWaveModes(hlm)

    def generate_td_waveform(self, **parameters):
        theta, phi = (
            parameters["inclination"],
            # parameters["phi_ref"]
            (np.pi / 2 - parameters["phi_ref"].value) * u.rad,
        )
        hlm = self.generate_td_modes(**parameters)
        hp, hc = hlm(theta, phi)
        hp, hc = TimeSeries(hp, name="hplus"), TimeSeries(hc, name="hcross")
        return hp, hc

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

        mode_array = parameters.pop("ModeArray",None)
        nr_data_file = parameters.pop('NumRelData')

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        f_start, dt = self.waveform_dict["f22_start"], self.waveform_dict["deltaT"]
        print('DBUG NR_hdf5 got dt=',dt)
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

        print('DBUG NR_hdf5 hp_lal', hp_lal, dir(hp_lal), dir(hp_lal.data))
        # print('DBUG NR_hdf5 hp_lal', hp_lal, dir(hp_lal), hp_lal.data.dt, hp_lal.data.t0)
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

class NRHybSur3dq8_CCE_gwsurr(NRHybSur3dq8_gwsurr):
    """
    Implements a toy wrapper for NRHybSur3dq8_CCE in the gwsurrogate package
    #TODO: Implement tapering
    """

    def __init__(self, **kwargs):
        self.sur = gwsurr.LoadSurrogate("NRHybSur3dq8_CCE")
        self._update_domains()
