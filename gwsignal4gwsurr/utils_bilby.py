"""
This file contains:

1. A ledger of surrogate models ready to use with bilby
2. Source functions used to generate fd waveforms
3. A parameter conversion function
4. A WaveformGenerator inherited class that routes code properly
"""

from dataclasses import dataclass, field
from bilby.gw.waveform_generator import WaveformGenerator
from bilby.gw.source import lal_binary_black_hole
from bilby.gw.conversion import (
    chirp_mass_and_mass_ratio_to_component_masses,
    convert_to_lal_binary_black_hole_parameters,
)
from .gwsurr import (
    NRHybSur3dq8_gwsurr,
    NRSur7dq4_gwsurr,
)
from .gwsurr_playground import (
    NRSur3dq8_Lev2_varenya_gwsurr,
    NRSur3dq8_Lev3_varenya_gwsurr,
    NRSur7dq4_LALSim_gwsurr,
)
import numpy as np
from typing import Optional, Any, Callable
from bilby.gw.conversion import bilby_to_lalsimulation_spins
from bilby.core.utils.constants import solar_mass
from astropy import units as u


# ================== ledger to hold surrogate model information ==================
@dataclass
class SurrogateWaveformProperties:
    precessing: bool
    wrapper: Callable
    _instance: Optional[Any] = field(default=None, init=False)
    eccentric: bool = False
    marginalization: bool = False

    # lazy load surrogate model instance
    @property
    def instance(self):
        if self._instance is None:
            self._instance = self.wrapper()
        return self._instance


surrogate_models = {
    "NRHybSur3dq8": SurrogateWaveformProperties(
        precessing=False, wrapper=NRHybSur3dq8_gwsurr
    ),
    "NRSur7dq4": SurrogateWaveformProperties(precessing=True, wrapper=NRSur7dq4_gwsurr),
    "NRSur3dq8_Lev2_varenya": SurrogateWaveformProperties(
        precessing=False, wrapper=NRSur3dq8_Lev2_varenya_gwsurr,marginalization=True
    ),
    "NRSur3dq8_Lev3_varenya": SurrogateWaveformProperties(
        precessing=False, wrapper=NRSur3dq8_Lev3_varenya_gwsurr,marginalization=True
    ),
    "NRSur7dq4_LALSim": SurrogateWaveformProperties(
        precessing=True, wrapper=NRSur7dq4_LALSim_gwsurr
    ),
}


# ===================== aligned spin frequency domain source model =====================
def gwsurrogate_binary_black_hole_aligned(
    freqs,
    mass1,
    mass2,
    spin1z,
    spin2z,
    distance,
    theta_jn,
    phi_ref,
    **waveform_arguments,
):
    f22_start = waveform_arguments["minimum_frequency"]

    approximant = waveform_arguments["waveform_approximant"]
    model = surrogate_models[approximant]
    gen = model.instance


    parameters = dict(
        mass1=mass1 * u.Msun,
        mass2=mass2 * u.Msun,
        spin1z=spin1z * u.dimensionless_unscaled,
        spin2z=spin2z * u.dimensionless_unscaled,
        distance=distance * u.Mpc,
        inclination=theta_jn * u.rad,
        phi_ref=(phi_ref) * u.rad,
        f22_start=f22_start * u.Hz,
        f22_ref=waveform_arguments["reference-frequency"] * u.Hz,
        f_max=waveform_arguments["maximum_frequency"] * u.Hz,
        deltaF=(freqs[1] - freqs[0]) * u.Hz,
    )

    # extra arguments to hand over to the model for marginalization
    extra_args = {}
    if model.marginalization:
        extra_args['noisy'] = waveform_arguments.get("noisy", False)
        extra_args['noise_level'] = waveform_arguments.get("noise_level", 1e-4)
        parameters['extra_args'] = extra_args
    try:
        hp_gwsignal, hc_gwsignal = gen.generate_fd_polarizations_from_td(
            **parameters
        )
    except Exception as e:
        if waveform_arguments["catch_waveform_errors"]:
            print(f"WARN surrogate wrapper failed to generate waveform: {e}")
            return None
        raise Exception(f"KILL surrogate wrapper failed to generate waveform: {e}")

    hp, hc = np.zeros_like(freqs, dtype=complex), np.zeros_like(freqs, dtype=complex)
    minimum_frequency = waveform_arguments["minimum_frequency"]
    maximum_frequency = waveform_arguments["maximum_frequency"]
    frequency_bounds = (freqs >= minimum_frequency) * (freqs <= maximum_frequency)

    if len(hp_gwsignal.data.data) > len(freqs):
        hp = hp_gwsignal.data.data[: len(hp)]
        hc = hc_gwsignal.data.data[: len(hc)]
    else:
        hp[: len(hp_gwsignal.data.data)] = hp_gwsignal.data.data
        hc[: len(hc_gwsignal.data.data)] = hc_gwsignal.data.data
    hp *= frequency_bounds
    hc *= frequency_bounds

    dt = 1 / hp_gwsignal.deltaF + (
        hp_gwsignal.epoch.gpsSeconds + hp_gwsignal.epoch.gpsNanoSeconds * 1e-9
    )
    time_shift = np.exp(-1j * 2 * np.pi * dt * freqs[frequency_bounds])
    hp[frequency_bounds] *= time_shift
    hc[frequency_bounds] *= time_shift

    return {"plus": hp, "cross": hc}


# ===================== precessing frequency domain source model =====================
def gwsurrogate_binary_black_hole_precessing(
    freqs,
    mass1,
    mass2,
    a_1,
    a_2,
    tilt_1,
    tilt_2,
    phi_12,
    phi_jl,
    distance,
    theta_jn,
    phi_ref,
    **waveform_arguments,
):
    f22_start = waveform_arguments["minimum_frequency"]

    # convert to cartesian spins
    # iota is angle between orbital angular momentum and line of sight (\theta_LN)
    iota, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = bilby_to_lalsimulation_spins(
        theta_jn=theta_jn,
        phi_jl=phi_jl,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        a_1=a_1,
        a_2=a_2,
        mass_1=mass1 * solar_mass,
        mass_2=mass2 * solar_mass,
        reference_frequency=waveform_arguments["reference-frequency"],
        phase=phi_ref,
    )

    approximant = waveform_arguments["waveform_approximant"]
    gen = surrogate_models[approximant].instance

    try:
        hp_gwsignal, hc_gwsignal = gen.generate_fd_polarizations_from_td(
            mass1=mass1 * u.Msun,
            mass2=mass2 * u.Msun,
            spin1x=spin1x * u.dimensionless_unscaled,
            spin1y=spin1y * u.dimensionless_unscaled,
            spin1z=spin1z * u.dimensionless_unscaled,
            spin2x=spin2x * u.dimensionless_unscaled,
            spin2y=spin2y * u.dimensionless_unscaled,
            spin2z=spin2z * u.dimensionless_unscaled,
            distance=distance * u.Mpc,
            inclination=iota * u.rad,
            phi_ref=phi_ref * u.rad,
            f22_start=f22_start * u.Hz,
            f22_ref=waveform_arguments["reference-frequency"] * u.Hz,
            f_max=waveform_arguments["maximum_frequency"] * u.Hz,
            deltaF=(freqs[1] - freqs[0]) * u.Hz,
        )
    except Exception as e:
        if waveform_arguments["catch_waveform_errors"]:
            print(f"PROG NRSur7dq4_wrapper failed to generate waveform: {e}")
            return None
        raise

    hp, hc = np.zeros_like(freqs, dtype=complex), np.zeros_like(freqs, dtype=complex)
    minimum_frequency = waveform_arguments["minimum_frequency"]
    maximum_frequency = waveform_arguments["maximum_frequency"]
    frequency_bounds = (freqs >= minimum_frequency) * (freqs <= maximum_frequency)

    if len(hp_gwsignal.data.data) > len(freqs):
        print(
            "WARN gwsurr waveform length exceeds requested frequency array length! Truncating waveform to fit frequency array."
        )
        hp = hp_gwsignal.data.data[: len(hp)]
        hc = hc_gwsignal.data.data[: len(hc)]
    else:
        hp[: len(hp_gwsignal.data.data)] = hp_gwsignal.data.data
        hc[: len(hc_gwsignal.data.data)] = hc_gwsignal.data.data
    hp *= frequency_bounds
    hc *= frequency_bounds

    dt = 1 / hp_gwsignal.deltaF + (
        hp_gwsignal.epoch.gpsSeconds + hp_gwsignal.epoch.gpsNanoSeconds * 1e-9
    )
    time_shift = np.exp(-1j * 2 * np.pi * dt * freqs[frequency_bounds])
    hp[frequency_bounds] *= time_shift
    hc[frequency_bounds] *= time_shift

    return {"plus": hp, "cross": hc}


# ===================== convert bilby params to gwsurrogate ones =====================
def parameter_conversion(parameters):
    # replace '-' with '_' in keys
    for k in parameters.keys():
        if "-" in k:
            k_ = k.replace("-", "_")
            parameters[k_] = parameters.pop(k)

    # initialize param dicts
    masses = {"mass1": None, "mass2": None}
    spins = {}
    extrinsic = {
        "distance": parameters["luminosity_distance"],
        "theta_jn": parameters["theta_jn"],
        "phi_ref": parameters["phase"],
    }

    # extract individual masses
    if "chirp_mass" in parameters.keys():
        masses["mass1"], masses["mass2"] = (
            chirp_mass_and_mass_ratio_to_component_masses(
                parameters["chirp_mass"], parameters["mass_ratio"]
            )
        )
    else:
        masses["mass1"] = parameters["mass_1"]
        masses["mass2"] = parameters["mass_2"]

    # spins
    if "chi_1" in parameters.keys():
        # aligned spin system specified with -1<chi_i<1
        spins["spin1z"] = parameters["chi_1"]
        spins["spin2z"] = parameters["chi_2"]

    elif "a_1" in parameters.keys():
        # precessing system; need to convert radial parameters to cartesian
        # since this conversion needs masses and f_ref, it is done in the wrapper function
        spins["a_1"] = parameters["a_1"]
        spins["a_2"] = parameters["a_2"]
        spins["tilt_1"] = parameters["tilt_1"]
        spins["tilt_2"] = parameters["tilt_2"]
        spins["phi_12"] = parameters["phi_12"]
        spins["phi_jl"] = parameters["phi_jl"]

    else:
        raise Exception("weird spin definition")

    converted_parameters = masses | spins | extrinsic

    keys = []
    for key in converted_parameters.keys():
        if key not in list(parameters.keys()):
            keys.append(key)

    return converted_parameters, keys


class SurrogateWaveformGenerator(WaveformGenerator):
    def __init__(self, **kwargs):
        approximant = kwargs["waveform_arguments"]["waveform_approximant"]

        model = surrogate_models.get(approximant, None)

        # if approximant isn't present in surrogate_models, revert to bilby defaults
        if model is None:
            print(f"INFO approximant {approximant} not implemented, reverting to bilby defaults")
            print(
                "INFO updating frequency_domain_source_model to", lal_binary_black_hole
            )
            print(
                "INFO updating parameter_conversion to",
                convert_to_lal_binary_black_hole_parameters,
            )

            kwargs["frequency_domain_source_model"] = lal_binary_black_hole
            kwargs["parameter_conversion"] = convert_to_lal_binary_black_hole_parameters

            super().__init__(**kwargs)
            return

        if model.precessing:
            kwargs["frequency_domain_source_model"] = (
                gwsurrogate_binary_black_hole_precessing
            )
        else:
            kwargs["frequency_domain_source_model"] = (
                gwsurrogate_binary_black_hole_aligned
            )
        kwargs["parameter_conversion"] = parameter_conversion

        super().__init__(**kwargs)
