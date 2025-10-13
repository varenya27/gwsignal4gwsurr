import numpy as np
from astropy import units as u
from .gwsurr import NRHybSur3dq8_gwsurr

gen = NRHybSur3dq8_gwsurr()
def NRHybSur3dq8_wrapper(freqs, mass1,mass2,spin1z,spin2z,distance,inclination,phi_ref,**waveform_arguments):
    if waveform_arguments['reference-frequency']<waveform_arguments['f-min']:
        print(f"DBUG fref {waveform_arguments['reference-frequency']} was lower than fmin {waveform_arguments['f-min']}! Setting fref=fmin")
        waveform_arguments['reference-frequency']=waveform_arguments['f-min']
    hp_gwsignal,hc_gwsignal =  gen.generate_fd_polarizations_from_td(
        mass1=mass1*u.Msun,
        mass2=mass2*u.Msun,
        spin1z=spin1z*u.dimensionless_unscaled,
        spin2z=spin2z*u.dimensionless_unscaled,
        distance=distance*u.Mpc,
        inclination=inclination*u.rad,
        phi_ref=(phi_ref)*u.rad,
        f22_start=waveform_arguments['f-min']*u.Hz,
        f22_ref=waveform_arguments['reference-frequency']*u.Hz,
        f_max = max(freqs)*u.Hz,
        deltaF=(freqs[1]-freqs[0])*u.Hz,
    )

    # VU: potential tc fix? cf bilby source.py#L647
    hp,hc = np.zeros_like(freqs,dtype=complex),np.zeros_like(freqs,dtype=complex)
    minimum_frequency = waveform_arguments['minimum_frequency']
    maximum_frequency = waveform_arguments['maximum_frequency']
    frequency_bounds = ((freqs >= minimum_frequency) * (freqs <= maximum_frequency))

    if len(hp_gwsignal.data.data)>len(freqs):
        hp = hp_gwsignal.data.data[:len(hp)]
        hc = hc_gwsignal.data.data[:len(hc)]
    else:
        hp[:len(hp_gwsignal.data.data)] = hp_gwsignal.data.data
        hc[:len(hc_gwsignal.data.data)] = hc_gwsignal.data.data
    hp *= frequency_bounds
    hc *= frequency_bounds

    dt = 1/hp_gwsignal.deltaF + (hp_gwsignal.epoch.gpsSeconds + hp_gwsignal.epoch.gpsNanoSeconds*1e-9)
    time_shift = np.exp(-1j * 2*np.pi * dt * freqs[frequency_bounds])
    hp[frequency_bounds] *= time_shift
    hc[frequency_bounds] *= time_shift

    return {'plus': hp, 'cross': hc}

