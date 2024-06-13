#custom libs
from utils.sclogging import log
import data.data_classes as dt
from utils.querypimms import *

#astronomical libraries
from astroquery.esa.xmm_newton import XMMNewton
from astropy.io import fits
from astropy.coordinates import Angle
from astropy import coordinates as coord
import astropy.units as u
import pyregion
import pyds9

#sas libraries
from pysas.wrapper import Wrapper as w

#scientific libraries
import numpy as np
import pandas as pd
import math
from scipy import stats
from scipy.optimize import curve_fit

#general libraries
import sys
import datetime
import os
import glob
import tarfile
from io import StringIO

#plotting libraries
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from screeninfo import get_monitors

def empty_folder(folder, delete=False):
    '''
    A function to empty a folder and optionally delete it afterwards.
    '''
    cwd = os.getcwd()
    os.chdir(folder)
    for file in glob.glob('*'):
        os.remove(file)
    os.chdir(cwd)
    if delete:
        os.rmdir(folder)
    return

def suppress_print_output(active=True, output_file='stdouterr.log'):
    '''
    A function to suppress print output to console and optionally write it to a file.
    '''
    if active:
        if output_file:
            sys.stdout = open(output_file, 'w')
            sys.stderr = open(output_file, 'w')
        else:
            sys.stdout = open(os.devnull, 'w')
            sys.stderr = open(os.devnull, 'w')
    else:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

def sas_command(cmd, args=[], dir=os.getcwd(), log_message=''):
    '''
    A function to run a SAS command with the given arguments in the specified directory.
    '''
    cwd = os.getcwd()
    os.chdir(dir)
    stdout_original = sys.stdout
    output_catturato = StringIO()
    in_error = False
    try:
        sys.stdout = output_catturato
        t = w(cmd, args)
        start = datetime.datetime.now()
        try:
            message = f'{log_message}' if log_message else ''
            log.info('Running command: ' + cmd + message)
            t.run()
            end = datetime.datetime.now()
            log.info(cmd + ' executed successfully in ' + str(end-start) + 's')
            log.debug('SAS command executed: %s' % cmd + ' inargs where ' + ', '.join(args))
        except Exception as e:
            in_error = True
            log.error('Error executing SAS command: %s' % cmd + ' inargs where ' + ', '.join(args))
    finally:
        sys.stdout = stdout_original
    output_catturato.seek(0)
    output = output_catturato.read()
    if in_error:
        log.warning('Error output: ' + output)
    else:
        log.debug('Output: ' + output)
    os.chdir(cwd)
    return output

def download_pps(observation, target_dir):
    '''
    A function to download the PPS data for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    if observation.specific_data.pps:
        log.debug('PPS data already downloaded')
        return
    obsid = observation.obsid
    log.info(f"Downloading pps data for {obsid}")
    inargs = ['odfid='+obsid, 'level=PPS']
    sas_command('startsas', inargs, target_dir)
    try:
        os.remove(target_dir+"/startsas.log")
    except:
        pass


def download_odf(observation, obs_dir):
    '''
    A function to download the ODF data for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of XMMNewtonObservation')
    if observation.specific_data.odf:
        log.debug('ODF data already downloaded')
        return
    cwd = os.getcwd()
    obsid = observation.obsid
    log.info(f"Downloading odf data for {obsid}")
    os.chdir(obs_dir)
    log.info(f"Downloading odf data for {obsid}")
    XMMNewton.download_data(obsid, level='ODF', extension='FTZ', extension_list=None, verbose=True)
    os.chdir(obs_dir)
    odf_dir = obs_dir+'odf/'
    os.makedirs(odf_dir, exist_ok=True)
    log.debug(f"Extracting odf data to {odf_dir}")
    
    #extract obsid.tar.gz to odf folder
    tar = tarfile.open(obsid+'.tar.gz')
    tar.extractall(odf_dir)
    tar.close()
    os.remove(obsid+'.tar.gz')

    #extract .tar files in odf folder
    os.chdir(odf_dir)
    for file in glob.glob('*.TAR'):
        print(file)
        tar = tarfile.open(file)
        tar.extractall(odf_dir)
        tar.close()
        os.remove(file)
    os.chdir(cwd)

def get_summary_file(observation):
    '''
    A function to get the summary file for a given observation with info about the instruments and exposures.
    '''
    cwd = os.getcwd()
    obs_dir = observation.dir_path
    pps_dir = obs_dir+'pps/'
    if not os.path.exists(pps_dir):
        raise FileNotFoundError('PPS directory not found')
    os.chdir(pps_dir)
    summary_file = glob.glob('P'+observation.obsid+'OBX000SUMMAR0000.HTM')[0]
    summary_file = pps_dir+'/'+summary_file
    os.chdir(cwd)
    return summary_file

def get_observation_info(observation):
    '''
    A function to get the observation info from the summary file and add it to the observation object.
    '''
    cwd = os.getcwd()
    obs_dir = observation.dir_path
    pps_dir = obs_dir+'pps/'
    if not os.path.exists(pps_dir):
        raise FileNotFoundError(f'PPS directory {pps_dir} not found')
    
    summary_file = get_summary_file(observation)
    log.debug('Summary file is: ' + summary_file)

    tables = pd.read_html(summary_file)
    obs_info = tables[3]
    inst_info = tables[4]
    exp_conf_info = tables[5]

    log.debug('Observation info: ' + str(obs_info))
    observation.start_date = obs_info[2][3]
    observation.end_date = obs_info[2][4]

    log.debug('Instrument info: ' + str(inst_info))
    instruments = {}
    for index, row in inst_info.iterrows():
        if row['Active'] == 'N':
            instruments[row['Instrument']] = {'active':False, 'modes':[], 'filters':[], 'exposure':0}
        elif row['Active'] == 'Y':
            instruments[row['Instrument']] = {'active':True, 'modes':[], 'filters':[], 'exposure':0}
    
    log.debug('Exposure configuration info: ' + str(exp_conf_info))
    for index, row in exp_conf_info.iterrows():
        if row['Inst.'] in instruments.keys() and instruments[row['Inst.']]['active'] and 'Diagnostic' not in str(row['Mode']) and row['Mode'] != 'Offset' and row['Mode'] != 'Diagnostic' and row['Mode'] != 'UNDEFINED':
            instruments[row['Inst.']]['modes'].append(str(row['Mode']))
            instruments[row['Inst.']]['filters'].append(str(row['Filter']))
            instruments[row['Inst.']]['exposure'] += row['Total Duration']
    
    for key in instruments.keys():
        instrument = dt.XMMNewtonInstrument(name=key, active=instruments[key]['active'], modes=list(set(instruments[key]['modes'])), filters=list(set(instruments[key]['filters'])), exposure=instruments[key]['exposure'])
        observation.specific_data.add_instrument(instrument)

    os.chdir(cwd)
    return

def cif_and_sumsas(observation, overwrite=False):
    '''
    A function to generate the CIF and SUM.SAS files for a given observation and set the environment variables.
    '''
    if observation.specific_data.cif and observation.specific_data.sumsas and not overwrite:
        os.environ['SAS_CCF'] = observation.specific_data.cif
        os.environ['SAS_ODF'] = observation.specific_data.sumsas
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        if observation.specific_data.odf:
            os.chdir(obs_dir)
            if not overwrite:
                try:
                    cif = glob.glob('ccf.cif')[0]
                    cif = obs_dir+cif
                    observation.specific_data.cif = cif
                except:
                    log.debug('CCF file not found, generating it')  
                    odf_dir = obs_dir+'odf/'
                    os.environ['SAS_ODF'] = odf_dir
                    try:
                        sas_command('cifbuild', [], obs_dir)
                    except:
                        log.error('Error while running cifbuild')
                        return
                    cif = glob.glob('ccf.cif')[0]
                    cif = obs_dir+cif
                    observation.specific_data.cif = cif
                try:
                    sumsas = glob.glob('*SUM.SAS')[0]
                    sumsas = obs_dir+sumsas
                    observation.specific_data.sumsas = sumsas
                except:
                    log.debug('SUM.SAS file not found, generating it')
                    os.environ['SAS_CCF'] = observation.specific_data.cif
                    try:
                        sas_command('odfingest', [], obs_dir)
                    except:
                        log.error('Error while running odfingest')
                        return
                    sumsas = glob.glob('*SUM.SAS')[0]
                    sumsas = obs_dir+sumsas
                    observation.specific_data.sumsas = sumsas
                os.chdir(cwd)
                return
            else:
                odf_dir = obs_dir+'odf/'
                os.environ['SAS_ODF'] = odf_dir
                try:
                    sas_command('cifbuild', [], obs_dir)
                except:
                    log.error('Error while running cifbuild')
                    return
                cif = glob.glob('ccf.cif')[0]
                cif = obs_dir+cif
                observation.specific_data.cif = cif
                os.environ['SAS_CCF'] = observation.specific_data.cif
                try:
                    sas_command('odfingest', [], obs_dir)
                except:
                    log.error('Error while running odfingest')
                    return
                sumsas = glob.glob('*SUM.SAS')[0]
                sumsas = obs_dir+sumsas
                observation.specific_data.sumsas = sumsas
                os.chdir(cwd)
        else:
            log.error('ODF directory not found, download it first')
            return

def get_observation_reduction_folders(observation):
    '''
    A function to get the reduction folders for a given observation and add them to the observation object.
    '''
    if observation.specific_data.instruments is None:
        try:
            get_observation_info(observation)
        except:
            log.error('Ob info absent and error while getting it')
            return
    obs_dir = observation.dir_path
    for instrument in observation.specific_data.active_instruments:
        if instrument.name == 'EMOS1' or instrument.name == 'EMOS2':
            emproc_dir = obs_dir+'emproc/'
            if os.path.exists(emproc_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('emproc', emproc_dir))
            emchain_dir = obs_dir+'emchain/'
            if os.path.exists(emchain_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('emchain', emchain_dir))
        if instrument.name == 'EPN':
            epnproc_dir = obs_dir+'epproc/'
            if os.path.exists(epnproc_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epproc', epnproc_dir))
            epnchain_dir = obs_dir+'epchain/'
            if os.path.exists(epnchain_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epchain', epnchain_dir))
            epnchainoot_dir = obs_dir+'epchainoot/'
            if os.path.exists(epnchainoot_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epchainoot', epnchain_dir))
        if instrument.name == 'RGS1' or instrument.name == 'RGS2':
            rgsproc_dir = obs_dir+'rgsproc/'
            if os.path.exists(rgsproc_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('rgsproc', rgsproc_dir))
        if instrument.name == 'OM':
            omichain_dir = obs_dir+'omichain/'
            if os.path.exists(omichain_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('omichain', omichain_dir))
            omfchain_dir = obs_dir+'omfchain/'
            if os.path.exists(omfchain_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('omfchain', omfchain_dir))
            omgchain_dir = obs_dir+'omgchain/'
            if os.path.exists(omgchain_dir):
                observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('omgchain', omgchain_dir))
    return

def emchain(observation, overwrite=False):
    '''
    A function to run the SAS emchain command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'EMOS1' not in inst_list and 'EMOS2' not in inst_list:
        log.debug('Neither EMOS1 nor EMOS2 active')
        return
    if 'emchain' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('emchain already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('emchain', [], tmp_dir, f' ({observation.obsid})')
        except:
            log.error(f'Error while executing emchain for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        emchain_dir = obs_dir+'emchain'
        os.makedirs(emchain_dir, exist_ok=True)
        empty_folder(emchain_dir)
        os.system(f'mv {tmp_dir}/* {emchain_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('emchain', emchain_dir))
        os.chdir(cwd)
    return

def emproc(observation, overwrite=False):
    '''
    A function to run the SAS emproc command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'EMOS1' not in inst_list and 'EMOS2' not in inst_list:
        log.debug('Neither EMOS1 nor EMOS2 active')
        return
    if 'emproc' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('emproc already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('emproc', [], tmp_dir, f' ({observation.obsid})')
        except:
            log.error(f'Error while executing emproc for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        emproc_dir = obs_dir+'emproc'
        os.makedirs(emproc_dir, exist_ok=True)
        empty_folder(emproc_dir)
        os.system(f'mv {tmp_dir}/* {emproc_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('emproc', emproc_dir))
        os.chdir(cwd)
    return

def epchain(observation, overwrite=False):
    '''
    A function to run the SAS epchain command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'EPN' not in inst_list:
        log.debug('EPN is not active')
        return
    if 'epchain' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug(f'epchain already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('epchain', ['runradmonfix=N'], tmp_dir, f' ({observation.obsid})')
        except:
            log.error(f'Error while executing epchain for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        epchain_dir = obs_dir+'epchain/'
        os.makedirs(epchain_dir, exist_ok=True)
        empty_folder(epchain_dir)
        os.system(f'mv {tmp_dir}/* {epchain_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epchain', epchain_dir))
        os.chdir(cwd)
    return

def epchainoot(observation, overwrite=False):
    '''
    A function to run the SAS epchainoot command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'EPN' not in inst_list:
        log.debug('EPN is not active')
        return
    if 'epchainoot' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug(f'epchainoot already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('epchain', ['runradmonfix=N', 'withoutoftime=true'], tmp_dir, f' ({observation.obsid})')
        except:
            log.error(f'Error while executing epchainoot for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        epchainoot_dir = obs_dir+'epchainoot/'
        os.makedirs(epchainoot_dir, exist_ok=True)
        empty_folder(epchainoot_dir)
        os.system(f'mv {tmp_dir}/* {epchainoot_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epchainoot', epchainoot_dir))
        os.chdir(cwd)
    return

def epproc(observation, overwrite=False):
    '''
    A function to run the SAS epproc command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'EPN' not in inst_list:
        log.debug('EPN is not active')
        return
    if 'epproc' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('epproc already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('epproc', [], tmp_dir, f'observation {observation.obsid}')
        except:
            log.error(f'Error while executing epproc for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        epproc_dir = obs_dir+'epproc/'
        os.makedirs(epproc_dir, exist_ok=True)
        empty_folder(epproc_dir)
        os.system(f'mv {tmp_dir}/* {epproc_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epproc', epproc_dir))
        os.chdir(cwd)
    return

def epprocoot(observation, overwrite=False):
    '''
    A function to run the SAS epprocoot command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'EPN' not in inst_list:
        log.debug('EPN is not active')
        return
    if 'epprocoot' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('epprocoot already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('epproc', ['withoutoftime=true'], tmp_dir, f' ({observation.obsid})')
        except:
            log.error(f'Error while executing epprocoot for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        epprocoot_dir = obs_dir+'epprocoot/'
        os.makedirs(epprocoot_dir, exist_ok=True)
        empty_folder(epprocoot_dir)
        os.system(f'mv {tmp_dir}/* {epprocoot_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('epprocoot', epprocoot_dir))
        os.chdir(cwd)
    return

def rgsproc(observation, overwrite=False):
    '''
    A function to run the SAS rgsproc command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'RGS1' not in inst_list and 'RGS2' not in inst_list:
        log.debug('Neither RGS1 nor RGS2 active')
        return
    if 'rgsproc' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('rgsproc already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('rgsproc', [], tmp_dir, f'observation {observation.obsid}')
        except:
            log.error(f'Error while executing rgsproc for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        rgsproc_dir = obs_dir+'rgsproc/'
        os.makedirs(rgsproc_dir, exist_ok=True)
        empty_folder(rgsproc_dir)
        os.system(f'mv {tmp_dir}/* {rgsproc_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('rgsproc', rgsproc_dir))
        os.chdir(cwd)
    return

def omichain(observation, overwrite=False):
    '''
    A function to run the SAS omichain command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'OM' not in inst_list:
        log.debug('OM is not active')
        return
    om_filters = observation.specific_data.get_instrument_by_name('OM').filters
    if not any(filter in om_filters for filter in ['U', 'UVW1', 'UVM2', 'UVW2', 'White', 'V', 'B']):
        log.debug('No valid filter found for omichain')
        return
    om_modes = observation.specific_data.get_instrument_by_name('OM').modes
    if not any(mode in om_modes for mode in ['Image']):
        log.debug('No image mode found')
        return
    if 'omichain' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('omichain already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('omichain', [], tmp_dir, f'observation {observation.obsid}')
        except:
            log.error(f'Error while executing omichain for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        omichain_dir = obs_dir+'omichain/'
        os.makedirs(omichain_dir, exist_ok=True)
        empty_folder(omichain_dir)
        os.system(f'mv {tmp_dir}/* {omichain_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('omichain', omichain_dir))
        os.chdir(cwd)
    return

def omfchain(observation, overwrite=False):
    '''
    A function to run the SAS omfchain command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'OM' not in inst_list:
        log.debug('OM is not active')
        return
    om_filters = observation.specific_data.get_instrument_by_name('OM').filters
    if not any(filter in om_filters for filter in ['U', 'UVW1', 'UVM2', 'UVW2', 'White', 'V', 'B']):
        log.debug('No valid filter found for omfchain')
        return
    om_modes = observation.specific_data.get_instrument_by_name('OM').modes
    if not any(mode in om_modes for mode in ['Fast']):
        log.debug('No fast mode found')
        return
    if 'omfchain' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('omfchain already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('omfchain', [], tmp_dir, f'observation {observation.obsid}')
        except:
            log.error(f'Error while executing omfchain for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        omfchain_dir = obs_dir+'omfchain/'
        os.makedirs(omfchain_dir, exist_ok=True)
        empty_folder(omfchain_dir)
        os.system(f'mv {tmp_dir}/* {omfchain_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('omfchain', omfchain_dir))
        os.chdir(cwd)
    return

def omgchain(observation, overwrite=False):
    '''
    A function to run the SAS omgchain command for a given observation.
    '''
    if not isinstance(observation, dt.Observation):
        raise TypeError('observation must be an instance of Observation')
    inst_list = [instrument.name for instrument in observation.specific_data.active_instruments]
    if 'OM' not in inst_list:
        log.debug('OM is not active')
        return
    om_filters = observation.specific_data.get_instrument_by_name('OM').filters
    if not any(filter in om_filters for filter in ['Grism 1', 'Grism 2']):
        log.debug('No valid filter found for omgchain')
        return
    if 'omgchain' in [reduction_folder.name for reduction_folder in observation.specific_data.reduction_folders] and not overwrite:
        log.debug('omgchain already executed')
        return
    else:
        cwd = os.getcwd()
        obs_dir = observation.dir_path
        tmp_dir = obs_dir+'tmp/'
        os.makedirs(tmp_dir, exist_ok=True)
        empty_folder(tmp_dir)
        try:
            cif_and_sumsas(observation)
            sas_command('omgchain', [], tmp_dir, f'observation {observation.obsid}')
        except:
            log.error(f'Error while executing omgchain for observation {observation.obsid}')
            empty_folder(tmp_dir)
            return
        omgchain_dir = obs_dir+'omgchain/'
        os.makedirs(omgchain_dir, exist_ok=True)
        empty_folder(omgchain_dir)
        os.system(f'mv {tmp_dir}/* {omgchain_dir}')
        empty_folder(tmp_dir, delete=True)
        observation.specific_data.add_reduction_folder(dt.XMMNewtonReduction('omgchain', omgchain_dir))
        os.chdir(cwd)
    return

def get_emos_evli_from_reduction(observation):
    '''
    A function to get the EMOS event list info for a given observation.
    '''
    avaiable_red = observation.specific_data.get_reduction_folders_names()
    if not any(red in avaiable_red for red in ['emchain', 'emproc']):
        log.debug('EMOS product are not available or not reduced') 
        return
    cwd = os.getcwd()
    obs_dir = observation.dir_path
    avaiable_emos_red = [red for red in avaiable_red if 'emchain' in red or 'emproc' in red]
    summary_file = get_summary_file(observation)
    tables = pd.read_html(summary_file)
    exp_conf_info = tables[5]
    if not os.path.exists(obs_dir+'images/'):
        os.makedirs(obs_dir+'images/')
    images_dir = obs_dir+'images/'
    for red in avaiable_emos_red:
        if red == 'emchain':
            red_dir = obs_dir+red+'/'
            os.chdir(red_dir)
            for file in glob.glob('*M1*MIEVLI*.FIT'):
                sensor = 'EMOS1'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('M1')[1].split('MIEVLI')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                image = dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration)
                observation.specific_data.add_image(image)
            for file in glob.glob('*M2*MIEVLI*.FIT'):
                sensor = 'EMOS2'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('M2')[1].split('MIEVLI')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            # for file in glob.glob('*'):
            #     if file not in [image.file for image in observation.specific_data.images]:
            #         os.remove(file)
        if red == 'emproc':
            red_dir = obs_dir+red+'/'
            os.chdir(red_dir)
            for file in glob.glob('*EMOS1_*_ImagingEvts.ds'):
                sensor = 'EMOS1'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('EMOS1_')[1].split('_ImagingEvts.ds')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            for file in glob.glob('*EMOS2_*_ImagingEvts.ds'):
                sensor = 'EMOS2'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('EMOS2_')[1].split('_ImagingEvts.ds')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            # for file in glob.glob('*'):
            #     if file not in [image.file for image in observation.specific_data.images]:
            #         os.remove(file)
    os.chdir(cwd)
            
def get_epn_evli_from_reduction(observation):
    '''
    A function to get the EPN event list info for a given observation.
    '''
    avaiable_red = observation.specific_data.get_reduction_folders_names()
    if not any(red in avaiable_red for red in ['epchain', 'epproc']):
        log.debug('EPN product are not available or not reduced') 
        return
    cwd = os.getcwd()
    obs_dir = observation.dir_path
    avaiable_epn_red = [red for red in avaiable_red if 'epchain' in red or 'epproc' in red]
    summary_file = get_summary_file(observation)
    tables = pd.read_html(summary_file)
    exp_conf_info = tables[5]
    if not os.path.exists(obs_dir+'images/'):
        os.makedirs(obs_dir+'images/')
    images_dir = obs_dir+'images/'
    for red in avaiable_epn_red:
        if red == 'epchain':
            red_dir = obs_dir+red+'/'
            os.chdir(red_dir)
            for file in glob.glob('*PN*PIEVLI*.FIT'):
                sensor = 'EPN'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('PN')[1].split('PIEVLI')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            # for file in glob.glob('*'):
            #     if file not in [image.file for image in observation.specific_data.images]:
            #         os.remove(file)
        if red == 'epchainoot':
            red_dir = obs_dir+red+'/'
            os.chdir(red_dir)
            for file in glob.glob('*PN*OOEVLI*.FIT'):
                sensor = 'EPN'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('PN')[1].split('OOEVLI')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            for file in glob.glob('*'):
                if file not in [image.file for image in observation.specific_data.images]:
                    os.remove(file)
        if red == 'epproc':
            red_dir = obs_dir+red+'/'
            os.chdir(red_dir)
            for file in glob.glob('*EPN*ImagingEvts.ds'):
                sensor = 'EPN'
                file = f'{file}'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('EPN_')[1].split('_ImagingEvts.ds')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            # for file in glob.glob('*'):
            #     if file not in [image.file for image in observation.specific_data.images]:
            #         os.remove(file)
        if red == 'epprocoot':
            red_dir = obs_dir+red+'/'
            os.chdir(red_dir)
            for file in glob.glob('*EPN*ImagingEvts.ds'):
                sensor = 'EPN'
                file = f'{file}'
                #rename file into file+OOT
                os.rename(file, file.split('.ds')[0]+'_OOT.ds')
                file = file.split('.ds')[0]+'_OOT.ds'
                os.system(f'cp {file} {obs_dir}images/')
                file_path = images_dir+file
                exp_number = file.split('EPN_')[1].split('_ImagingEvts_OOT.ds')[0]
                row = exp_conf_info[(exp_conf_info['Inst.'] == sensor) & (exp_conf_info['Exp Id'] == exp_number)]
                mode = row['Mode'].values[0]
                filter = row['Filter'].values[0]
                duration = row['Total Duration'].values[0]
                observation.specific_data.add_image(dt.XMMNewtonImage(red, sensor, exp_number, file, file_path, mode, filter, duration))
            for file in glob.glob('*'):
                if file not in [image.file for image in observation.specific_data.images]:
                    os.remove(file)
    os.chdir(cwd)

#TODO: add rgsproc, omichain, omfchain, omgchain get_info functions

def open_ds9(evn_file):
    '''
    A function to open a ds9 window with the given event file.
    '''
    os.environ['LD_LIBRARY_PATH'] = ''
    os.environ['PATH'] = '/opt/SAS/xmmsas_20230412_1735/binextra:/opt/SAS/xmmsas_20230412_1735/bin:/opt/SAS/xmmsas_20230412_1735/bin/devel:/opt/HEAsoft/x86_64-pc-linux-gnu-libc2.35/bin:/home/andrea/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin'

    monitors = get_monitors()
    if len(monitors) > 1:  # If there are multiple monitors, open DS9 on the secondary screen
        monitor = monitors[1]  # Secondary screen
        d = pyds9.DS9(start=f"-geometry {int(monitor.height//(9/8))}x{int(monitor.height//(9/8))}+{int(monitor.x)+int((monitor.width-int(monitor.height//(9/8)))//2)}+{int((monitor.height-int(monitor.height//(9/8)))//4)}")
    else:
        monitor = monitors[0]  # Primary screen
        d = pyds9.DS9(start=f"-geometry {int(monitor.width//2)}x{int(monitor.width//2)}+{int(monitor.width//2)}+{int((monitor.height-monitor.width//2)//4)}")

    d.set(f'file {evn_file}')
    d.set('frame center')
    d.set('bin factor 64')
    d.set('scale log')
    d.set('scale limits 1.5 150')
    d.set('zoom to fit')
    return d

def open_ds9s(evn_files):
    '''
    A function to open a ds9 window with the given event files.
    '''
    os.environ['LD_LIBRARY_PATH'] = ''
    os.environ['PATH'] = '/opt/SAS/xmmsas_20230412_1735/binextra:/opt/SAS/xmmsas_20230412_1735/bin:/opt/SAS/xmmsas_20230412_1735/bin/devel:/opt/HEAsoft/x86_64-pc-linux-gnu-libc2.35/bin:/home/andrea/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin'

    monitors = get_monitors()
    if len(monitors) > 1:  # If there are multiple monitors, open DS9 on the secondary screen
        monitor = monitors[1]  # Secondary screen
        d = pyds9.DS9(start=f"-geometry {int(monitor.height//(9/8))}x{int(monitor.height//(9/8))}+{int(monitor.x)+int((monitor.width-int(monitor.height//(9/8)))//2)}+{int((monitor.height-int(monitor.height//(9/8)))//4)}")
    else:
        monitor = monitors[0]  # Primary screen
        d = pyds9.DS9(start=f"-geometry {int(monitor.width//2)}x{int(monitor.width//2)}+{int(monitor.width//2)}+{int((monitor.height-monitor.width//2)//4)}")
    
    for i, file in enumerate(evn_files):
        if i > 0:
            d.set('frame new')
        try:
            d.set(f'file {file}')
        except:
            log.error(f'Error while opening {file}')
            continue
        d.set('frame center')
        d.set('bin factor 64')
        d.set('scale zscale')
        d.set('zoom to fit')
    d.set('tile')
    return d

def copy_to(file_path, dest_dir):
    '''
    A function to copy a file to a destination directory.
    '''
    os.system(f'cp {file_path} {dest_dir}')

def copy_all_images_to_wd(observation):
    '''
    A function to copy all images to the working directory.
    '''
    obs_dir = observation.dir_path
    wd_dir = obs_dir+'wd/'
    os.makedirs(wd_dir, exist_ok=True)
    for image in observation.specific_data.images:
        copy_to(image.file_path, wd_dir)

def get_ds9_regions(d):
    '''
    A function to get the regions drawn in a ds9 window.
    '''
    d.set(f'regions save tmpregion.reg')
    try:
        regions = pyregion.open('tmpregion.reg')
        os.remove('tmpregion.reg')
    except:
        log.warning('No regions found or error in saving file or removing it')
        regions = []
    return regions

def description_ds9_regions(regions):
    '''
    A function to add a description to each region in a list of regions.
    '''
    regions_dict = {}
    for i in range(len(regions)):
        region_str = str(regions[i].name).upper() + '(' + ','.join([str(x) for x in regions[i].coord_list]) + ')'
        # description = input(f'Description for region {i}: ')
        # regions_dict[f'{i}'] = {'description':description, 'region':region_str}
        regions_dict[f'{i}'] = {'region':region_str}
    return regions_dict

def plotLC_flaringbkg(plt,threshold,fileName):
    '''
    A function to plot a light curve with a threshold.
    ''' 
    if fileName != "NOT FOUND":
        fitsFile = fits.open(fileName)   
        prihdu   = fitsFile[1].header
        if ('CUTVAL' in prihdu):
            threshold = prihdu['CUTVAL']

        cols    = fitsFile[1].columns

        colName = None
        for i,x  in enumerate(cols.names):
            if "RATE" in x:
                colName = cols.names[i]
            if "COUNTS" in x:
                colName = cols.names[i]        

        if colName is None:
            raise ValueError("Neither 'RATE' nor 'COUNTS' found in column names")

        data = fitsFile[1].data   
        
        xdata = data.field('TIME') - min(data.field('TIME')) # extract the x data column
        start_time = data.field('TIME')[0]
        start_time = datetime.datetime(1998, 1, 1) + datetime.timedelta(seconds=start_time)
        ydata = data.field(colName)

        xmax=np.amax(xdata) 
        xmin=np.amin(xdata) 

        plt.plot(xdata,ydata) # plot the data

        if colName == 'RATE':
            plt.title("Flaring particle background lightcurve")
            plt.xlabel("Time (s)")
            plt.ylabel("Cts/s")
        else:
            plt.title("Lightcurve")
            plt.xlabel("Time (s)")
            plt.ylabel("Counts")

        if (threshold != 'None'):
            if (colName == 'COUNTS'):
                threshold=float(threshold)*100.

            y2data = [threshold]*len(xdata)            
            plt.plot(xdata,y2data)
            plt.text(xmin+0.1*(xmax-xmin), threshold+0.01*threshold, str(threshold)+" cts/sec", ha='center')        

        fitsFile.close()   
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xlabel(f'(s) - start time: {start_time}')
    else:
        print("File not found "+fileName+"\n")

def find_proc_region_files(data_dir, obsid):
    '''
    A function to find the region files in the regions directory for a given observation.
    '''
    cwd = os.getcwd()
    obs_dir = data_dir+'/'+obsid
    regions_dir = obs_dir+'/regions'
    regions_dict = {}
    os.chdir(regions_dir)
    for file in glob.glob('*proc.reg'):
            regions_dict[file] = str(os.path.join(regions_dir, file))
    os.chdir(cwd)
    return regions_dict

def find_chain_region_files(data_dir, obsid):
    '''
    A function to find the region files in the regions directory for a given observation.
    '''
    cwd = os.getcwd()
    obs_dir = data_dir+'/'+obsid
    regions_dir = obs_dir+'/regions'
    regions_dict = {}
    os.chdir(regions_dir)
    for file in glob.glob('*chain.reg'):
            regions_dict[file] = str(os.path.join(regions_dir, file))
    os.chdir(cwd)
    return regions_dict

##check
# def select_sources_regions(obsid, data_dir, method='proc'):
#     obs_dir = data_dir+'/'+obsid
#     os.chdir(obs_dir)
#     if method == 'proc':
#         evn_files = find_proc_products(obsid, data_dir)
#     elif method == 'chain':
#         evn_files = find_chain_products(obsid, data_dir)
#     else:
#         log.error('Invalid method')
    
#     wd_dir = obs_dir+'/wd'
#     os.mkdirs(wd_dir, exist_ok=True)
#     copy_evn_files_to(data_dir, obsid, evn_files, wd_dir)

#     os.chdir(wd_dir)
#     d1 = open_ds9s(evn_files)

#     #user will draw a region for each frame, then ask him to confirm
#     log.info('Asking user to draw sources regions')
#     input('Draw a region for each frame to select sources, then press enter to continue')

#     #check if the regions directory exists, if not create it
#     os.mkdirs('../regions', exist_ok=True)

#     #save all regions as ds9 regions files
#     for i in range(len(evn_files)):
#         j=i+1
#         d1.set(f'frame {j}')
#         d1.set(f'regions save source_region_{list(evn_files.keys())[i]}.reg')
#         os.system(f'mv source_region_{list(evn_files.keys())[i]}.reg ../regions')
#         log.info(f"Region saved as {os.path.join(os.getcwd(), f'source_region_{list(evn_files.keys())[i]}.reg')}")
#     os.chdir(obs_dir)
#     d1.set('exit')

#     log.info('Sources region files saved in working directory')

# ##check
# def select_background_regions(obsid, data_dir, method='proc'):
#     obs_dir = data_dir+'/'+obsid
#     os.chdir(obs_dir)
#     if method == 'proc':
#         evn_files = find_proc_products(obsid, data_dir)
#     elif method == 'chain':
#         evn_files = find_chain_products(obsid, data_dir)
#     else:
#         log.error('Invalid method')
    
#     wd_dir = obs_dir+'/wd'
#     os.mkdirs(wd_dir, exist_ok=True)
#     copy_evn_files_to(data_dir, obsid, evn_files, wd_dir)

#     os.chdir(wd_dir)
#     d = open_ds9s(evn_files)

#     #user will draw a region for each frame, then ask him to confirm
#     log.info('Asking user to draw background regions')
#     input('Draw a region for each frame to select background, then press enter to continue')

#     #save all regions as ds9 regions files
#     for i in range(len(evn_files)):
#         j=i+1
#         d.set(f'frame {j}')
#         d.set(f'regions save background_region_{list(evn_files.keys())[i]}.reg')
#         os.system(f'mv background_region_{list(evn_files.keys())[i]}.reg ../regions')
#         log.info(f"Region saved as {os.path.join(os.getcwd(), f'background_region_{list(evn_files.keys())[i]}.reg')}")
#     d.set('exit')
#     log.info('Background region files saved in working directory')   

def filter_by_expression(input_evt_table, output_evt_table, expression, dir):
    '''
    A function to filter an event table by an expression.
    '''
    inargs = []
    inargs.append('table='+input_evt_table)
    inargs.append('withfilteredset=yes')
    inargs.append('filteredset='+output_evt_table)
    inargs.append(f'expression={expression}')
    log.debug('Filtering by expression: %s' % expression)
    sas_command('evselect', inargs, dir, f'filtering by expression {expression}')

def filter_by_region(input_evt_table, region, add_str=''):
    '''
    A function to filter an event table by a region.
    '''
    if not isinstance(region,pyregion.Shape):
        raise TypeError('region must be an instance of pyregion.Shape')
    regione = str(region.name).upper() + '(' + ','.join([str(x) for x in region.coord_list]) + ')'
    region_str = f'(X,Y) IN {regione}'
    log.debug('Filtering by region: %s' % region_str)
    region_str_name = str(region.name).upper() + '_' + '_'.join([str(x) for x in region.coord_list])
    output_evt_table = input_evt_table.replace('.ds', f'_{region_str_name}.ds')
    if add_str != '' and add_str not in output_evt_table:
        output_evt_table = output_evt_table.replace('.ds', f'_{add_str}.ds')
    filter_by_expression(input_evt_table, output_evt_table, region_str)
    return output_evt_table
        
def filter_by_pi(input_evt_table, pi_range, add_str=''):
    '''
    A function to filter an event table by PI range.
    '''
    expression = f'PI>={pi_range[0]} && PI<={pi_range[1]}'
    log.debug('Filtering by PI: PI IN %s' % str(pi_range))
    output_evt_table = input_evt_table.replace('.ds', f'_PI{pi_range[0]}_{pi_range[1]}.ds')
    if add_str != '' and add_str not in output_evt_table:
        output_evt_table = output_evt_table.replace('.ds', f'_{add_str}.ds')
    filter_by_expression(input_evt_table, output_evt_table, expression)
    return output_evt_table

def filter_by_pattern_and_flag(input_evt_table, pattern, flag, add_str=''):
    '''
    A function to filter an event table by pattern and flag.
    '''
    expression = f'PATTERN<={pattern} && FLAG=={flag}'
    log.debug('Filtering by pattern and flag: PATTERN<=%s && FLAG==%s' % (pattern, flag))
    output_evt_table = input_evt_table.replace('.ds', f'_P{pattern}_F{flag}.ds')
    if add_str != '' and add_str not in output_evt_table:
        output_evt_table = output_evt_table.replace('.ds', f'_{add_str}.ds')
    filter_by_expression(input_evt_table, output_evt_table, expression)
    return output_evt_table

def filter_by_rate(input_evt_table, rate, add_str=''):
    '''
    A function to filter an event table by rate.
    '''
    expression = f'RATE<{rate}'
    log.debug('Filtering by rate: RATE<%s' % rate)
    output_evt_table = input_evt_table.replace('.ds', f'_R{rate}.ds')
    if add_str != '' and add_str not in output_evt_table:
        output_evt_table = output_evt_table.replace('.ds', f'_{add_str}.ds')
    filter_by_expression(input_evt_table, output_evt_table, expression)
    return output_evt_table

##TODO
# def filter_by_gti(input_evt_table, output_evt_table, input_gti):

# def tabgtigen(input_lc, output_gti, expression):

# def tabgtigen_for_bkg_filtering(input_lc, output_gti, rate):

def light_curve(input_evt_table, output_lc, time_bin, dir, expression=None, filt_set_name=None):
    '''
    A function to create a light curve from an event table.
    '''
    inargs = []
    inargs.append('table='+input_evt_table)
    inargs.append('withrateset=yes')
    inargs.append('rateset='+output_lc)
    inargs.append('maketimecolumn=yes')
    inargs.append('makeratecolumn=yes')
    inargs.append('timebinsize='+str(time_bin))
    if expression:
        inargs.append(f'expression={expression}')
    if filt_set_name:
        inargs.append('withfilteredset=yes')
        inargs.append('filteredset='+filt_set_name)
    log.debug(f'Creating light curve for file {input_evt_table} with time bin {time_bin}')
    sas_command('evselect', inargs, dir, f'light curve for {input_evt_table}')
    return

def epiclccorr(input_evt_table, input_lc, output_lc, dir, bkg_lc=None):
    '''
    A function to correct a light curve.
    '''
    inargs = []
    inargs.append('srctslist='+input_lc)
    inargs.append('eventlist='+input_evt_table)
    inargs.append('applyabsolutecorrections=yes')
    if bkg_lc:
        inargs.append('withbkgset=yes')
        inargs.append('bkgtslist='+bkg_lc)
    inargs.append('outset='+output_lc)
    log.debug(f'Running epiclccorr for {input_lc}')
    sas_command('epiclccorr', inargs, dir)
    return output_lc

def get_max_rate(input_lc):
    '''
    A function to get the maximum rate from a light curve.
    '''
    light_curve = fits.open(input_lc)
    rate = light_curve[1].data['RATE']
    max_rate = max(rate)
    return max_rate

def draw_lc(input_lc, label, save=True, ylims=None, scale=None, color=None, gti=None):
    '''
    A function to draw a light curve.
    '''
    light_curve = fits.open(input_lc)
    time = light_curve[1].data['TIME']
    rate = light_curve[1].data['RATE']
    if scale:
        rate = rate*scale
    error = light_curve[1].data['ERROR']
    if scale:
        error = error*scale
    start_time = light_curve[1].data['TIME'][0]
    start_time = datetime.datetime(1998, 1, 1) + datetime.timedelta(seconds=start_time)
    time0 = time[0]
    time = time - time0
    plt.errorbar(time, rate, yerr=error, markersize=1, elinewidth=0.5, capsize=1, capthick=0.5, color=color, label=label)
    plt.xlabel(f'(s) - start time: {start_time}')
    plt.ylabel('(cts/s)')
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    if ylims:
        plt.ylim(ylims[0], ylims[1])
    else:  
        max_rate = max(rate)
        ylims = [0, max_rate+0.1*max_rate]
        log.info(f'Max rate: {max_rate}')
        plt.ylim(ylims[0], ylims[1])
    
    if label:
        plt.legend()

    if gti:
        #gti is a list of time intervals, if gti is not None, change the colour of the points inside the gtis
        gti = fits.open(gti)
        gti = gti[1].data
        for i in range(len(gti)):
            plt.axvspan(gti[i][0]-time0, gti[i][1]-time0, alpha=0.2, color='grey')

    if save:
        plt.savefig(input_lc+'.png')

    return time0, ylims

def draw_lc_corr(input_lc, label, save=True, ylims=None, scale=None, color=None, gti=None):
    '''
    A function to draw a corrected light curve.
    '''
    light_curve = fits.open(input_lc)
    time = light_curve[1].data['TIME']
    rate = light_curve[1].data['RATE']
    if scale:
        rate = rate*scale
    error = light_curve[1].data['ERROR']
    if scale:
        error = error*scale
    start_time = light_curve[1].data['TIME'][0]
    start_time = datetime.datetime(1998, 1, 1) + datetime.timedelta(seconds=start_time)
    time0 = time[0]
    time = time - time0
    ax1 = plt.gca()
    ax1.errorbar(time, rate, yerr=error, markersize=1, elinewidth=0.5, capsize=1, capthick=0.5, color=color, label=label)
    ax1.spines['left'].set_position(('axes', 0))
    ax1.set_ylabel('Y', color='black')
    plt.xlabel(f'(s) - start time: {start_time}')
    plt.ylabel('(cts/s)')
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    if ylims:
        plt.ylim(ylims[0], ylims[1])
    else:  
        max_rate = max(rate[~np.isnan(rate)])
        ylims = [0, max_rate+0.1*max_rate]
        plt.ylim(ylims[0], ylims[1])
    
    if label:
        plt.legend()

    if gti:
        #gti is a list of time intervals, if gti is not None, change the colour of the points inside the gtis
        gti = fits.open(gti)
        gti = gti[1].data
        for i in range(len(gti)):
            plt.axvspan(gti[i][0]-time0, gti[i][1]-time0, alpha=0.2, color='grey')

    if save:
        plt.savefig(input_lc+'.png')

    return ax1

def draw_hr(input_lc, color=None):
    '''
    A function to draw a hardness ratio.
    '''
    light_curve = fits.open(input_lc)
    time = light_curve[1].data['TIME']
    rate = light_curve[1].data['RATE']
    error = light_curve[1].data['ERROR']
    #remove negative values
    error[error<0] = 0
    start_time = light_curve[1].data['TIME'][0]
    start_time = datetime.datetime(1998, 1, 1) + datetime.timedelta(seconds=start_time)
    time = time - time[0]
    #get the current axis
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    if not color:
        color = 'black'
    ax2.errorbar(time, rate, yerr=error, markersize=5, elinewidth=0.5, capsize=1, capthick=0.5, linestyle='-', marker='o', color=color)
    ax2.spines['right'].set_position(('axes', 1))
    ax2.set_ylabel('Y', color='black')
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.ylabel('Hardness ratio')
    plt.ylim(-1, 1)

def draw_temp(input_temp, color=None):
    '''
    A function to draw a temperature.
    '''
    light_curve = fits.open(input_temp)
    time = light_curve[1].data['TIME']
    temperature = light_curve[1].data['RATE']
    error = light_curve[1].data['ERROR']
    #remove negative values
    error[error<0] = 0
    start_time = light_curve[1].data['TIME'][0]
    start_time = datetime.datetime(1998, 1, 1) + datetime.timedelta(seconds=start_time)
    time = time - time[0]
    #get the current axis
    ax1 = plt.gca()
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 20))
    ax3.set_ylabel('Y3 Label', color='violet')
    if not color:
        color = 'black'
    ax3.errorbar(time, temperature, yerr=error, markersize=5, elinewidth=0.5, capsize=1, capthick=0.5, linestyle='-', marker='o', color=color)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.ylabel('Temperature')
    plt.ylim(6, 9)

def hardness_ratio(infile1, infile2, outfile):
    '''
    A function to calculate the hardness ratio between two light curves.
    '''
    lc1 = fits.open(infile1)
    lc2 = fits.open(infile2)
    rate1 = lc1[1].data['RATE']
    rate2 = lc2[1].data['RATE']
    error1 = lc1[1].data['ERROR']
    error2 = lc2[1].data['ERROR']
    try:
        hr = (rate1-rate2)/(rate1+rate2)
        hr_error = ((error1**2+error2**2)/(rate1+rate2)**2)**0.5
    except:
        log.error('Error while calculating hardness ratio')
        return
    lc1[1].data['RATE'] = hr
    lc1[1].data['ERROR'] = hr_error
    lc1.writeto(outfile, overwrite=True)
    return outfile

def dsplot(input):
    '''
    A function to plot a table.
    '''
    os.system(f'dsplot table={input}')

def preliminary_filter_mos(input_evli, output_evli, dir, energy_range=[300, 1000]):
    '''
    Preliminary filter for MOS event lists.
    '''
    inargs = []
    inargs.append('table='+input_evli)
    inargs.append('withimageset=yes')
    inargs.append('imageset='+output_evli)
    inargs.append('filtertype=expression')
    inargs.append(f'expression="(PI in [{energy_range[0]}:{energy_range[1]}])&&(PATTERN<=12)&&((FLAG & 0x766aa000)==0)"')
    inargs.append('ignorelegallimits=yes')
    inargs.append('imagebinning=imageSize')
    inargs.append('xcolumn=DETX')
    inargs.append('ximagesize=780')
    inargs.append('ximagemax=19500')
    inargs.append('ximagemin=-19499')
    inargs.append('ycolumn=DETY')
    inargs.append('yimagesize=780')
    inargs.append('yimagemax=19500')
    inargs.append('yimagemin=-19499')
    sas_command('evselect', inargs, dir, ' (preliminary_filter_mos)')

def preliminary_filter_pn(input_evli, output_evli, dir, energy_range=[300, 1000]):
    '''
    Preliminary filter for PN event lists.
    '''
    inargs = []
    inargs.append('table='+input_evli)
    inargs.append('withimageset=yes')
    inargs.append('imageset='+output_evli)
    inargs.append('filtertype=expression')
    inargs.append(f'expression="(PI in [{energy_range[0]}:{energy_range[1]}])&&(PATTERN <= 4)&&(#XMMEA_EP)"')
    inargs.append('ignorelegallimits=yes')
    inargs.append('imagebinning=imageSize')
    inargs.append('xcolumn=DETX')
    inargs.append('ximagesize=780')
    inargs.append('ximagemax=19500')
    inargs.append('ximagemin=-19499')
    inargs.append('ycolumn=DETY')
    inargs.append('yimagesize=780')
    inargs.append('yimagemax=19500')
    inargs.append('yimagemin=-19499')
    sas_command('evselect', inargs, dir, ' (preliminary_filter_pn)')

def check_anomalous_ccds(image_path, dir):
    '''
    A function to check for anomalous CCDs.
    '''
    inargs = []
    inargs.append('eventfile='+image_path)
    inargs.append('keepcorner=no')
    sas_command('emanom', inargs, dir)

def generate_gti(input_evt_table, output_gti, threshold, dir):
    '''
    A function to generate a GTI from an event table based on a rate threshold.
    '''
    inargs = []
    inargs.append('table='+input_evt_table)
    inargs.append(f'expression="(RATE<{threshold})"')
    inargs.append('gtiset='+output_gti)
    sas_command('tabgtigen', inargs, dir)

def merge_gtis(gti1, gti2, output_gti):
    '''
    A function to merge two GTIs.
    '''
    copy_to(gti1, output_gti)
    #open gti files with astropy.io fits
    gti1 = fits.open(gti1)
    gti2 = fits.open(gti2)
    gti1_data = gti1[1].data
    gti2_data = gti2[1].data
    #gtis are lists of time intervals, i want to intersect gti1 with gti2
    #gti1 is the reference gti
    new_gti = []
    for i in range(len(gti1_data)):
        for j in range(len(gti2_data)):
            if gti2_data[j][0] > gti1_data[i][0] and gti2_data[j][1] < gti1_data[i][1]:
                new_gti.append([gti2_data[j][0], gti2_data[j][1]])
    out_gti = fits.open(output_gti, mode='update')
    new_gti_array = np.array(new_gti, dtype=[('START', 'f8'), ('STOP', 'f8')])
    out_gti[1].data = new_gti_array
    out_gti.flush()
    out_gti.close()
    return output_gti

def draw_lc_diff(lc_1, lc_2, color='red'):
    '''
    A function to draw the difference between two light curves.
    '''
    lc_1 = fits.open(lc_1)
    lc_2 = fits.open(lc_2)
    time_1 = lc_1[1].data['TIME']
    rate_1 = lc_1[1].data['RATE']
    error_1 = lc_1[1].data['ERROR']
    time_2 = lc_2[1].data['TIME']
    rate_2 = lc_2[1].data['RATE']
    error_2 = lc_2[1].data['ERROR']
    time_diff = time_1
    rate_diff = rate_1 - rate_2
    start_time = lc_1[1].data['TIME'][0]
    start_time = datetime.datetime(1998, 1, 1) + datetime.timedelta(seconds=start_time)
    time0 = time_diff[0]
    time_diff = time_diff - time0
    #draw the light curve in blue witouth error bars
    plt.plot(time_diff, rate_diff, color=color)
    #close the files
    lc_1.close()
    lc_2.close()

def eregionanalyse(imageset, srcexp, bkgexp, dir):
    '''
    A function to perform eregionanalyse and get the optimal source region radius.
    '''
    inargs = []
    inargs.append('imageset='+imageset)
    inargs.append('srcexp='+srcexp)
    inargs.append('backexp='+bkgexp)
    sas_command('eregionanalyse', inargs, dir)

def spectrum(in_set, out_set=None, out_spectrum=None, expression=None, binsize=None, min=None, max=None, dir=None):
    '''
    A function to create a spectrum from an event table. It saves the filtered data and the spectrum data.
    '''
    inargs = []
    inargs.append('table='+in_set)
    inargs.append('energycolumn=PI')
    if out_set is not None:
        inargs.append('withfilteredset=yes')
        inargs.append('filteredset='+out_set)
    else:
        inargs.append('withfilteredset=no')
    inargs.append('keepfilteroutput=yes')
    inargs.append('filtertype=expression')
    inargs.append(f'expression={expression}')
    inargs.append('withspectrumset=yes')
    inargs.append('spectrumset='+out_spectrum)
    inargs.append('spectralbinsize='+str(binsize))
    inargs.append('withspecranges=yes')
    inargs.append('specchannelmin='+str(min))
    inargs.append('specchannelmax='+str(max))
    # # inargs.append('updateexposure=yes')
    # inargs.append('filterexposure=yes')
    sas_command('evselect', inargs, dir, ' (spectrum)')

def backscale(in_spectrum, in_set, dir, ignoreoutoffov='yes', badpixelresolution=0.5):
    '''
    A function to backscale a spectrum and get the ratio between the source and background regions.
    '''
    inargs = []
    inargs.append('spectrumset='+in_spectrum)
    inargs.append('badpixlocation='+in_set)
    inargs.append(f'ignoreoutoffov={ignoreoutoffov}')
    inargs.append(f'badpixelresolution={badpixelresolution}')
    sas_command('backscale', inargs, dir)

def rmfgen(in_spectrum, out_rmf, dir):
    '''
    A function to create an RMF file from a spectrum.
    '''
    inargs = []
    inargs.append('spectrumset='+in_spectrum)
    inargs.append('rmfset='+out_rmf)
    inargs.append('withenergybins=yes')
    inargs.append('threshold=0')
    inargs.append('energymin=0.0')
    inargs.append('energymax=12.0')
    inargs.append('nenergybins=100')
    #inargs.append('correctforpileup=yes')
    sas_command('rmfgen', inargs, dir)

def arfgen(in_spectrum, in_set, in_rmf, out_arf, dir, computebackscales=False, withsourceposition=False, X=None, Y=None, sourcecoords='pos'):
    '''
    A function to create an ARF file from a spectrum.
    '''
    inargs = []
    inargs.append('spectrumset='+in_spectrum)
    inargs.append('arfset='+out_arf)
    inargs.append('withrmfset=yes')
    inargs.append('rmfset='+in_rmf)
    inargs.append('badpixlocation='+in_set)
    inargs.append('detmaptype=psf')
    inargs.append('eegridfactor=2')
    inargs.append('extendedsource=no')
    inargs.append('modelee=yes')
    #inargs.append('applyxcaladjustment=yes')
    inargs.append('setbackscale=yes')
    if withsourceposition:
        inargs.append('withsourcepos=yes')
        inargs.append(f'sourcecoords="{sourcecoords}"')
        inargs.append(f'sourcex={X}')
        inargs.append(f'sourcey={Y}')
    sas_command('arfgen', inargs, dir)

def specgroup(in_spectrum, in_rmf, in_arf, in_background, out_group, mincounts, oversample, dir):
    '''
    A function to group a spectrum.
    '''
    inargs = []
    inargs.append('spectrumset='+in_spectrum)
    inargs.append('rmfset='+in_rmf)
    inargs.append('arfset='+in_arf)
    if in_background is not None:
        inargs.append('backgndset='+in_background)
    inargs.append('groupedset='+out_group)
    inargs.append('mincounts='+str(mincounts))
    inargs.append('oversample='+str(oversample))
    sas_command('specgroup', inargs, dir)

def get_backscale_value(in_spectrum):
    '''
    A function to get the backscale value from a spectrum.
    '''
    spectrum = fits.open(in_spectrum)
    backscale = spectrum[1].header['BACKSCAL']
    return backscale

def set_backscale_value(in_spectrum, value):
    '''
    A function to set the backscale value in a spectrum.
    '''
    spectrum = fits.open(in_spectrum, mode='update')
    spectrum[1].header['BACKSCAL'] = value
    spectrum.flush()
    spectrum.close()

def get_average_error(lc_file):
    '''
    A function to get the average error from a light curve.
    '''
    lc = fits.open(lc_file)
    error = lc[1].data['ERROR']
    average_error = np.average(error)
    return average_error

def get_median_error(lc_file):
    '''
    A function to get the median error from a light curve.
    '''
    lc = fits.open(lc_file)
    error = lc[1].data['ERROR']
    average_error = np.median(error)
    return average_error

def merge_gti(gti1, gti2, mode="intersection"):
    '''
    A function to merge two GTIs.
    '''
    file1 = fits.open(gti1)[1]
    file2 = fits.open(gti2)[1]
    start1 = file1.header['TSTART']
    stop1 = file1.header['TSTOP']
    start2 = file2.header['TSTART']
    stop2 = file2.header['TSTOP']
    log.debug(f'GTI1: {start1} - {stop1}')
    log.debug(f'GTI2: {start2} - {stop2}')
    intervals1 = [x for x in file1.data]
    log.debug(f'Intervals1: {intervals1}')
    intervals2 = [x for x in file2.data]
    log.debug(f'Intervals2: {intervals2}')
    if mode == "union":
        if start1 < start2:
            start = start1
        else:
            start = start2
        if stop1 > stop2:
            stop = stop1
        else:
            stop = stop2
        intervals = []
        for i in range(len(intervals1)):
            intervals.append(intervals1[i])
        for i in range(len(intervals2)):
            intervals.append(intervals2[i])
        intervals = sorted(intervals, key=lambda x: x[0])
        new_intervals = []
        for i in range(len(intervals)):
            if i == 0:
                new_intervals.append(intervals[i])
            else:
                if intervals[i][0] <= new_intervals[-1][1]:
                    if intervals[i][1] > new_intervals[-1][1]:
                        new_intervals[-1] = (new_intervals[-1][0], intervals[i][1])
                else:
                    new_intervals.append(intervals[i])
        intervals = new_intervals
    elif mode == "intersection":
        if start1 < start2:
            start = start2
        else:
            start = start1
        if stop1 > stop2:
            stop = stop2
        else:
            stop = stop1
        intervals = []
        for i in range(len(intervals1)):
            for j in range(len(intervals2)):
                if intervals1[i][0] < intervals2[j][1] and intervals1[i][1] > intervals2[j][0]:
                    if intervals1[i][0] < intervals2[j][0]:
                        start = intervals2[j][0]
                    else:
                        start = intervals1[i][0]
                    if intervals1[i][1] > intervals2[j][1]:
                        stop = intervals2[j][1]
                    else:
                        stop = intervals1[i][1]
                    intervals.append((start, stop))
    else:
        log.error('Invalid mode')
        sys.exit()
    log.debug(f'Intervals: {intervals}')
    return start, stop, intervals

def invert_gti(start, stop, intervals):
    '''
    A function to invert a GTI.
    '''
    log.debug(f'Inverting GTI: {start} - {stop}')
    log.debug(f'Intervals: {intervals}')
    new_intervals = []
    if intervals[0][0] > start:
        new_intervals.append((start, intervals[0][0]))
    for i in range(len(intervals)):
        if i == 0:
            new_intervals.append((intervals[i][1], intervals[i+1][0]))
        elif i == len(intervals) - 1:
            new_intervals.append((intervals[i][1], stop))
        else:
            new_intervals.append((intervals[i][1], intervals[i+1][0]))
    log.debug(f'New intervals: {new_intervals}')
    return start, stop, new_intervals

def save_gti(start, stop, intervals, file):
    '''
    A function to save a GTI to a file.
    '''
    #save intervals in a simple text file
    with open(file, 'w') as f:
        f.write(str(start) + ' ' + str(stop) + '\n')
        for interval in intervals:
            f.write(str(interval[0]) + ' ' + str(interval[1]) + '\n')
    log.info(f'GTI saved in {file}')

def load_gti(file):
    '''
    A function to load a GTI from a file.
    '''
    #load intervals from a simple text file
    with open(file, 'r') as f:
        lines = f.readlines()
    start, stop = lines[0].strip().split()
    start = float(start)
    stop = float(stop)
    intervals = []
    for line in lines[1:]:
        interval = line.strip().split()
        intervals.append((float(interval[0]), float(interval[1])))
    return start, stop, intervals

def draw_intervals(intervals, start):
    '''
    A function to draw intervals on a plot.
    '''
    intervals = [(float(x[0]) - float(start), float(x[1]) - float(start)) for x in intervals]
    for i in range(len(intervals)):
        plt.axvspan(intervals[i][0], intervals[i][1], color='gray', alpha=0.3)

def draw_spectrum(file, range=None, softband=None, hardband=None, scale=None):
    '''
    A function to draw a spectrum.
    '''
    spectrum = fits.open(file)
    binsize = spectrum[1].header['SPECDELT']
    energy = spectrum[1].data['CHANNEL']*binsize/1000
    counts = spectrum[1].data['COUNTS']
    if scale:
        counts = counts*scale
    plt.plot(energy, counts)
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    if range:
        plt.xlim(range[0], range[1])
    if softband:
        plt.axvspan(softband[0], softband[1], color='green', alpha=0.3)
    if hardband:
        plt.axvspan(hardband[0], hardband[1], color='red', alpha=0.3)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

def get_first_time(input_lc):
    '''
    A function to get the first time from a light curve.
    '''
    light_curve = fits.open(input_lc)
    time = light_curve[1].data['TIME']
    time0 = time[0]
    return time0

def get_last_time(input_lc):
    '''
    A function to get the last time from a light curve.
    '''
    light_curve = fits.open(input_lc)
    time = light_curve[1].data['TIME']
    lasttime = time[-1]
    return lasttime

def epatplot(evli, dir, bgset=None):
    '''
    A function to plot the EPAT plot.
    '''
    inargs = []
    inargs.append('set='+evli)
    name_is = evli.replace('.FIT','')
    inargs.append(f'plotfile={name_is}-epatplot.ps')
    if bgset:
        inargs.append('backgroundset='+bgset)
    output = sas_command('epatplot', inargs, dir)
    return output

def draw_counts_histogram(input_lc, bins=25):
    '''
    A function to draw a histogram of counts.
    '''
    if isinstance(input_lc, str):
        light_curve = fits.open(input_lc)
        counts = light_curve[1].data['RATE']
    if isinstance(input_lc, list):
        light_curve1 = fits.open(input_lc[0])
        light_curve2 = fits.open(input_lc[1])
        counts1 = light_curve1[1].data['RATE']
        counts2 = light_curve2[1].data['RATE']
        counts = np.concatenate((counts1, counts2))
    #remove nan values
    counts = counts[~np.isnan(counts)]
    plt.hist(counts, bins=bins)
    plt.xlabel('Counts')
    plt.ylabel('Frequency')

def save_region_in_header(file, region, name, description=None):
    '''
    A function to save a region in the header of a file.
    '''
    #open the file
    hdul = fits.open(file, mode='update')
    #get the header
    header = hdul[1].header
    if description:
        region = region+'/'+description
    #add the region to the header
    header[name] = region
    #close the file
    hdul.close()

def remove_region_in_header(file, region, name):
    '''
    A function to remove a region in the header of a file.
    '''
    #open the file
    hdul = fits.open(file, mode='update')
    #get the header
    header = hdul[1].header
    #remove the region from the header
    del header[name]
    #close the file
    hdul.close()

def get_region_from_header(file, name, description=False):
    '''
    A function to get a region from the header of a file.
    '''
    #open the file
    hdul = fits.open(file)
    #get the header
    header = hdul[1].header
    #get the region from the header
    try:
        value = header[name]
    except:
        log.error(f'No region found in header with name {name}')
        return ''
    if description:
        description = value.split('/')[1]
        region = value.split('/')[0]
    else:
        region = value
    #close the file
    hdul.close()
    if description:
        return region, description
    return region

def add_excised_radius(file, linen, sensor, value):
    #aggiungi il valore in fondo alla riga linen del file .txt
    #il carattere di separazione è \t
    with open(file, 'r') as f:
        lines = f.readlines()
    line = lines[linen]
    line = line.strip() + '\t' + str(sensor) + '\t' + str(value) + '\n'
    lines[linen] = line
    with open(file, 'w') as f:
        f.writelines(lines)
    return
    
    






###################
# XSPEC FUNCTIONS #
###################

