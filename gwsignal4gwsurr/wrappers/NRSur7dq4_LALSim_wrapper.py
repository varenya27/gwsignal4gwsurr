from .gwsurr import NRSur7dq4_gwsurr
import numpy as np
from bilby.gw.conversion import bilby_to_lalsimulation_spins
from bilby.core.utils.constants import solar_mass

gen = NRSur7dq4_gwsurr()


def NRSur7dq4_wrapper(
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

    # if waveform_arguments['reference-frequency']<waveform_arguments['f-min']:
    #     print(f"DBUG fref {waveform_arguments['reference-frequency']} was lower than fmin {waveform_arguments['f-min']}! Setting fref=fmin")
    #     waveform_arguments['reference-frequency']=waveform_arguments['f-min']
    #
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
    # print('DBUG converted spins', iota, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z)

    # check if the discrepancies come at the strain generation stage
    import lalsimulation as lalsim
    import lal

    try:
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
            iota,
            phi_ref,
            0.0,
            0.0,
            0.0,
            freqs[1] - freqs[0],
            waveform_arguments["minimum_frequency"],
            waveform_arguments["maximum_frequency"],
            waveform_arguments["reference-frequency"],  # f_low, f_max, f_ref
            lal.CreateDict(),
            lalsim.GetApproximantFromString("NRSur7dq4"),
        )
    except Exception as e:
        print(f"PROG NRSur7dq4_wrapper failed to generate waveform: {e}")
        return None
    # print('DBUG gwsig lalsim used params', mass1,mass2,spin1x,spin1y,spin1z,spin2x,spin2y,spin2z,distance,iota,phi_ref,waveform_arguments['minimum_frequency'],max(freqs),waveform_arguments['reference-frequency'],freqs[1]-freqs[0])

    # print('DBUG epoch for gwsig lalsim =',hp_gwsignal.epoch.gpsSeconds , hp_gwsignal.epoch.gpsNanoSeconds)
    # return {'plus': hp_gwsignal.data.data, 'cross': hc_gwsignal.data.data}

    hp, hc = np.zeros_like(freqs, dtype=complex), np.zeros_like(freqs, dtype=complex)
    minimum_frequency = waveform_arguments["minimum_frequency"]
    maximum_frequency = waveform_arguments["maximum_frequency"]
    frequency_bounds = (freqs >= minimum_frequency) * (freqs <= maximum_frequency)
    # print('DBUG frequency bounds', frequency_bounds)

    if len(hp_gwsignal.data.data) > len(freqs):
        print(
            "WARN gwsurr waveform length exceeds requested frequency array length! Truncating waveform to fit frequency array."
        )
        hp = hp_gwsignal.data.data[: len(hp)]
        hc = hc_gwsignal.data.data[: len(hc)]
    else:
        # print('DBUG lalsim waveform lengths hp, hc=', len(hp_gwsignal.data.data), len(hc_gwsignal.data.data))
        hp[: len(hp_gwsignal.data.data)] = hp_gwsignal.data.data
        hc[: len(hc_gwsignal.data.data)] = hc_gwsignal.data.data
    hp *= frequency_bounds
    hc *= frequency_bounds

    dt = 1 / hp_gwsignal.deltaF + (
        hp_gwsignal.epoch.gpsSeconds + hp_gwsignal.epoch.gpsNanoSeconds * 1e-9
    )
    # print('DBUG gwsig lal hp_gwsignal.deltaF,epoch=', hp_gwsignal.deltaF, hp_gwsignal.epoch.gpsSeconds, hp_gwsignal.epoch.gpsNanoSeconds)
    # print('DBUG gwsig lal time shift dt=', dt)
    time_shift = np.exp(-1j * 2 * np.pi * dt * freqs[frequency_bounds])
    hp[frequency_bounds] *= time_shift
    hc[frequency_bounds] *= time_shift

    return {"plus": hp, "cross": hc}
