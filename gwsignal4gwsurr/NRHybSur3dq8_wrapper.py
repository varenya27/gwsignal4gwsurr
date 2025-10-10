from astropy import units as u
from .gwsurr import NRHybSur3dq8_gwsurr

gen = NRHybSur3dq8_gwsurr()
def NRHybSur3dq8_wrapper(freqs, mass1,mass2,spin1z,spin2z,distance,inclination,phi_ref,**waveform_arguments):
    if waveform_arguments['reference-frequency']<waveform_arguments['f-min']:
        print(f"DBUG fref {waveform_arguments['reference-frequency']} was lower than fmin {waveform_arguments['f-min']}! Setting fref=fmin")
        waveform_arguments['reference-frequency']=waveform_arguments['f-min']
    hp,hc =  gen.generate_fd_polarizations_from_td(
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
    return {'plus': hp.data, 'cross': hc.data}

