'''
templates for parameter conversion and waveform generator
use these two in any new wrapper and send to bilby
'''
from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_component_masses, convert_to_lal_binary_black_hole_parameters
from bilby.gw.waveform_generator import WaveformGenerator
from bilby.gw.source import lal_binary_black_hole
import numpy as np 

def parameter_conversion_aligned(parameters):
    # replace '-' with '_' in keys
    for k in parameters.keys():
        if '-' in k:
            k_ = k.replace('-','_')
            parameters[k_] = parameters.pop(k)

    masses = {'mass1':None, 'mass2':None}
    spins = {}
    extrinsic = {
        'distance':parameters['luminosity_distance'],
        'theta_jn':parameters['theta_jn'],
        'phi_ref':parameters['phase']
    }

    # extract individual masses
    if 'chirp_mass' in parameters.keys():
        masses['mass1'],masses['mass2'] = chirp_mass_and_mass_ratio_to_component_masses(parameters['chirp_mass'],parameters['mass_ratio'])
    else:
        masses['mass1'] = parameters['mass_1']
        masses['mass2'] = parameters['mass_2']

    # spins
    if 'chi_1' in parameters.keys():
        # print('DBUG converting spins')
        # aligned spin system specified with -1<chi_i<1
        spins['spin1z'] = parameters['chi_1']
        spins['spin2z'] = parameters['chi_2']

    elif 'a_1' in parameters.keys():
        # aligned spin system specified with 0<a_i<1 and tilt=[0,pi]
        spins['spin1z'] = parameters.pop('a_1')*np.cos(parameters.pop('tilt_1'))
        spins['spin2z'] = parameters.pop('a_2')*np.cos(parameters.pop('tilt_2'))
    
    else:
        raise Exception('weird spin definition')

    converted_parameters = masses | spins | extrinsic

    keys=[]
    for key in converted_parameters.keys():
        if key not in list(parameters.keys()):
            keys.append(key)
    return converted_parameters, keys

def parameter_conversion_precessing(parameters):
    # replace '-' with '_' in keys
    for k in parameters.keys():
        if '-' in k:
            k_ = k.replace('-','_')
            parameters[k_] = parameters.pop(k)

    masses = {'mass1':None, 'mass2':None}
    spins = {}
    extrinsic = {
        'distance':parameters['luminosity_distance'],
        'theta_jn':parameters['theta_jn'],
        'phi_ref':parameters['phase']
    }

    # extract individual masses
    if 'chirp_mass' in parameters.keys():
        masses['mass1'],masses['mass2'] = chirp_mass_and_mass_ratio_to_component_masses(parameters['chirp_mass'],parameters['mass_ratio'])
    else:
        masses['mass1'] = parameters['mass_1']
        masses['mass2'] = parameters['mass_2']

    # precessing system; need to convert radial parameters to cartesian
    # since this conversion needs masses and f_ref, it is done in the wrapper function
    spins['a_1']    = parameters['a_1']
    spins['a_2']    = parameters['a_2']
    spins['tilt_1'] = parameters['tilt_1']
    spins['tilt_2'] = parameters['tilt_2']
    spins['phi_12'] = parameters['phi_12']
    spins['phi_jl'] = parameters['phi_jl']

    converted_parameters = masses | spins | extrinsic

    keys=[]
    for key in converted_parameters.keys():
        if key not in list(parameters.keys()):
            keys.append(key)
    return converted_parameters, keys

parameter_conversions = [
    parameter_conversion_aligned,
    parameter_conversion_precessing
]

def get_waveform_generator(**kwargs):

    # for injections, revert to bilby defaults
    if kwargs['waveform_arguments']['waveform_approximant'] == 'NR_hdf5':
        print('INFO reverting to bilby defaults for injection')
        print('INFO updating frequency_domain_source_model to', lal_binary_black_hole)
        print('INFO updating parameter_conversion to', convert_to_lal_binary_black_hole_parameters)

        kwargs['frequency_domain_source_model'] = lal_binary_black_hole
        kwargs['parameter_conversion'] = convert_to_lal_binary_black_hole_parameters
        return WaveformGenerator(**kwargs)

    # bilby sometimes defaults to the inbuilt BBH parameter conversion function, which we don't want
    # [TODO] this is VERY hacky, need to figure out a better way to do this 
    if not kwargs['parameter_conversion'] in parameter_conversions:
        if '3d' in str(kwargs['frequency_domain_source_model']):
            print(f"PROG Updating parameter conversion function from {kwargs['parameter_conversion']} to {parameter_conversion_aligned}")
            kwargs['parameter_conversion'] = parameter_conversion_aligned
        elif '7d' in str(kwargs['frequency_domain_source_model']):
            print(f"PROG Updating parameter conversion function from {kwargs['parameter_conversion']} to {parameter_conversion_precessing}")
            kwargs['parameter_conversion'] = parameter_conversion_precessing
        else:
            raise Exception('Not sure what conversion function to assign for this waveform model')


    return WaveformGenerator(**kwargs)
