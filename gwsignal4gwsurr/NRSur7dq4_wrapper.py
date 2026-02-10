'''
    NOT IMPLEMENTED YET | DOES NOT WORK!
'''
from astropy import units as u
from .gwsurr import NRSur7dq4_gwsurr
import numpy as np
from bilby.gw.conversion import bilby_to_lalsimulation_spins
from bilby.core.utils.constants import solar_mass

gen = NRSur7dq4_gwsurr()
def NRSur7dq4_wrapper(freqs, mass1,mass2,a_1,a_2,tilt_1,tilt_2,phi_12, phi_jl,distance,theta_jn,phi_ref,**waveform_arguments):

    if 'f-min' in waveform_arguments.keys():
        f22_start = waveform_arguments['f-min']
    else:
        f22_start = waveform_arguments['minimum_frequency']

    # print('DBUG f22_start=', f22_start)
    if waveform_arguments['reference-frequency']<f22_start:
        print(f"DBUG fref {waveform_arguments['reference-frequency']} was lower than fmin {f22_start}! Setting fref=fmin")
        waveform_arguments['reference-frequency']=f22_start


    # convert to cartesian spins
    # iota is angle between orbital angular momentum and line of sight (\theta_LN)
    # print('DBUG converting', theta_jn, phi_jl, tilt_1, tilt_2, phi_12, a_1, a_2, mass1, mass2, waveform_arguments['reference-frequency'], phi_ref)
    iota, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = bilby_to_lalsimulation_spins(
        theta_jn = theta_jn,
        phi_jl =  phi_jl,
        tilt_1 = tilt_1,
        tilt_2 = tilt_2,
        phi_12 = phi_12,
        a_1 = a_1,
        a_2 = a_2,
        mass_1 = mass1*solar_mass,
        mass_2 = mass2*solar_mass,
        reference_frequency = waveform_arguments['reference-frequency'],
        phase = phi_ref
    )
    # print('DBUG got', iota, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z)

    if waveform_arguments['catch_waveform_errors']:
        try:
            hp_gwsignal,hc_gwsignal =  gen.generate_fd_polarizations_from_td(
                mass1=mass1*u.Msun,
                mass2=mass2*u.Msun,
                spin1x=spin1x*u.dimensionless_unscaled,
                spin1y=spin1y*u.dimensionless_unscaled,
                spin1z=spin1z*u.dimensionless_unscaled,
                spin2x=spin2x*u.dimensionless_unscaled,
                spin2y=spin2y*u.dimensionless_unscaled,
                spin2z=spin2z*u.dimensionless_unscaled,
                distance=distance*u.Mpc,
                inclination=iota*u.rad,
                phi_ref=phi_ref*u.rad,
                f22_start=f22_start*u.Hz,
                f22_ref=waveform_arguments['reference-frequency']*u.Hz,
                f_max = waveform_arguments['maximum_frequency']*u.Hz,
                deltaF=(freqs[1]-freqs[0])*u.Hz,
            )
        except Exception as e:
            print(f"PROG NRSur7dq4_wrapper failed to generate waveform: {e}")
            return None
    else:
        hp_gwsignal,hc_gwsignal =  gen.generate_fd_polarizations_from_td(
            mass1=mass1*u.Msun,
            mass2=mass2*u.Msun,
            spin1x=spin1x*u.dimensionless_unscaled,
            spin1y=spin1y*u.dimensionless_unscaled,
            spin1z=spin1z*u.dimensionless_unscaled,
            spin2x=spin2x*u.dimensionless_unscaled,
            spin2y=spin2y*u.dimensionless_unscaled,
            spin2z=spin2z*u.dimensionless_unscaled,
            distance=distance*u.Mpc,
            inclination=iota*u.rad,
            phi_ref=phi_ref*u.rad,
            f22_start=f22_start*u.Hz,
            f22_ref=waveform_arguments['reference-frequency']*u.Hz,
            f_max = waveform_arguments['maximum_frequency']*u.Hz,
            deltaF=(freqs[1]-freqs[0])*u.Hz,
        )

    # print('DBUG gwsig gwsurr used params', mass1,mass2,spin1x,spin1y,spin1z,spin2x,spin2y,spin2z,distance,iota,phi_ref,waveform_arguments['minimum_frequency'],max(freqs),waveform_arguments['reference-frequency'],freqs[1]-freqs[0])

    # print('DBUG epoch for gwsig gwsurr =',hp_gwsignal.epoch.gpsSeconds, hp_gwsignal.epoch.gpsNanoSeconds)
    # return {'plus': hp_gwsignal.data.data, 'cross': hc_gwsignal.data.data}

    hp,hc = np.zeros_like(freqs,dtype=complex),np.zeros_like(freqs,dtype=complex)
    minimum_frequency = waveform_arguments['minimum_frequency']
    maximum_frequency = waveform_arguments['maximum_frequency']
    frequency_bounds = ((freqs >= minimum_frequency) * (freqs <= maximum_frequency))
    # print('DBUG frequency bounds', frequency_bounds)

    if len(hp_gwsignal.data.data)>len(freqs):
        print('WARN gwsurr waveform length exceeds requested frequency array length! Truncating waveform to fit frequency array.')
        hp = hp_gwsignal.data.data[:len(hp)]
        hc = hc_gwsignal.data.data[:len(hc)]
    else:
        # print('DBUG gwsurr waveform lengths hp, hc=', len(hp_gwsignal.data.data), len(hc_gwsignal.data.data))
        hp[:len(hp_gwsignal.data.data)] = hp_gwsignal.data.data
        hc[:len(hc_gwsignal.data.data)] = hc_gwsignal.data.data
    hp *= frequency_bounds
    hc *= frequency_bounds

    dt = 1/hp_gwsignal.deltaF + (hp_gwsignal.epoch.gpsSeconds + hp_gwsignal.epoch.gpsNanoSeconds*1e-9)
    time_shift = np.exp(-1j * 2*np.pi * dt * freqs[frequency_bounds])
    hp[frequency_bounds] *= time_shift
    hc[frequency_bounds] *= time_shift

    return {'plus': hp, 'cross': hc}

