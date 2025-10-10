'''
    NOT IMPLEMENTED YET | DOES NOT WORK!
'''
from astropy import units as u
from .gwsurr import NRSur7dq4_gwsurr

from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_component_masses
from bilby.gw.waveform_generator import WaveformGenerator

gen = NRSur7dq4_gwsurr()
def NRSur7dq4_wrapper(freqs, mass1,mass2,spin1z,spin2z,distance,inclination,phi_ref,**waveform_arguments):
    if waveform_arguments['reference_frequency']<waveform_arguments['f_min']:
        print(f"DBUG fref {waveform_arguments['reference_frequency']} was lower than fmin {waveform_arguments['f_min']}! Setting fref=fmin")
        waveform_arguments['reference_frequency']=waveform_arguments['f_min']
    hp,hc =  gen.generate_fd_polarizations_from_td(
        mass1=mass1*u.Msun,
        mass2=mass2*u.Msun,
        spin1z=spin1z*u.dimensionless_unscaled,
        spin2z=spin2z*u.dimensionless_unscaled,
        distance=distance*u.Mpc,
        inclination=inclination*u.rad,
        phi_ref=(phi_ref)*u.rad,
        f22_start=waveform_arguments['f_min']*u.Hz,
        f22_ref=waveform_arguments['reference_frequency']*u.Hz,
        f_max = max(freqs)*u.Hz,
        deltaF=(freqs[1]-freqs[0])*u.Hz,
    )
    return {'plus': hp.data, 'cross': hc.data}

