'''
templates for parameter conversion and waveform generator
use these two in any new wrapper and send to bilby
'''
from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_component_masses
from bilby.gw.waveform_generator import WaveformGenerator

def parameter_conversion(parameters):
    mass_1,mass_2 = chirp_mass_and_mass_ratio_to_component_masses(parameters['chirp_mass'],parameters['mass_ratio'])
    converted_parameters = {
        'mass1':mass_1,
        'mass2':mass_2,
        'spin1z':parameters['chi_1'],
        'spin2z':parameters['chi_2'],
        'distance':parameters['luminosity_distance'],
        'inclination':parameters['theta_jn'],
        'phi_ref':parameters['phase'],
    }
    keys=[]
    for key in converted_parameters.keys():
        if key not in list(parameters.keys()):
            keys.append(key)
    return converted_parameters, keys

def get_waveform_generator(**kwargs):
    # bilby sometimes defaults to the inbuilt BBH parameter conversion function, which we don't want
    if not kwargs['parameter_conversion'] is parameter_conversion:
        print(f"PROG Updating parameter conversion function from {kwargs['parameter_conversion']} to {parameter_conversion}")
        kwargs['parameter_conversion']=parameter_conversion

    return WaveformGenerator(**kwargs)


