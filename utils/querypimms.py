#custom libs
from utils.sclogging import log

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
import datetime
import time
import os
import requests
import glob
import tarfile
import subprocess

#plotting libraries
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

def call_to_pimms(commands):
	'''
	Function to call PIMMS with a list of commands and return the output.
	'''
	process = subprocess.Popen("/opt/pimms4_13a/pimms", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	output, _ = process.communicate(input=b"".join(commands))
	process.stdin.close()
	process.stdout.close()
	process.stderr.close()
	return output.decode('utf-8')

# Inputs passed to command-line tool:
# FROM XMM PN THIN .3-1.0
# MODEL plasma 6.00 logT 0.8 17
# PLASMA MEKAL
# INSTRUMENT XMM PN THIN 1.0-2.0
# GO 1

def hr_from_pimms(inst, filt, nh=17, abundance=0.8, model="rs"):
	'''
	Function to calculate the hardness ratio from PIMMS.
	'''
	if inst not in ["pn", "mos"]:
		raise ValueError("Instrument must be 'pn' or 'mos'")
	if filt not in ["thin", "medium", "thick"]:
		raise ValueError("Filter must be 'thin', 'medium' or 'thick'")
	if model not in ["rs", "mekal"]:
		raise ValueError("Model must be 'rs' or 'mekal'")
	#from 5.6 to 8.5 in steps of 0.01
	logT_values = [
		round(5.6 + 0.05 * i, 2) for i in range(0, 65)
	]
	hr_values = []
	for value in logT_values:
		commands = [
			f"from xmm {inst} {filt} .3-1.0\n",
			f"model plasma {value} logT {abundance} {nh}\n",
			f"plasma {model}\n",
			f"instrument xmm {inst} {filt} 1.0-2.0\n",
			f"go 1\n",
		]
		#convert all commands to bytes
		commands = [command.encode('utf-8') for command in commands]
		output = call_to_pimms(commands)
		counts = float(output.split('* PIMMS predicts ')[1].split(' cps')[0])
		hr_values.append((counts*0.7-1)/(counts*0.7+1))
	hr_values = list(zip(hr_values, logT_values))
	return hr_values

def temperature_from_hr(hr_values, hr_target):
	'''
	Function to calculate the temperature from the hardness ratio.
	'''
	closest_temp = None
	min_diff = float('inf')
	for hr, temp in hr_values:
		diff = abs(hr - hr_target)
		if diff < min_diff:
			closest_temp = temp
			min_diff = diff

	return closest_temp  # Returning only the temperature part of the tuple

def convert_hr_to_temperature(hr_lc, inst, filt, nh=17, abundance=0.8, model="rs"):
	'''
	Function to convert a hardness ratio lightcurve to a temperature lightcurve.
	'''
	hr_values = hr_from_pimms(inst, filt, nh, abundance, model)
	#open file and get the hardness ratio
	hdu = fits.open(hr_lc)
	hr = hdu[1].data['RATE']
	hrerr = hdu[1].data['ERROR']
	#get the temperature
	temperature = []
	for h in hr:
		temperature.append(temperature_from_hr(hr_values, h))
	temperature = np.array(temperature)
	temperatureerr = np.zeros(len(temperature))
	for i in range(len(temperature)):
		try:
			min_temp = temperature_from_hr(hr_values, h - hrerr[i])
			max_temp = temperature_from_hr(hr_values, h + hrerr[i])
			temperatureerr[i] = (max_temp - min_temp) / 2
		except:
			temperatureerr[i] = 100
			log.warning(f"Error in temperature error calculation at index {i}")
	#save a new file with the temperature
	hdu[1].data['RATE'] = temperature
	hdu[1].data['ERROR'] = temperatureerr
	hdu.writeto(hr_lc[:-5] + "temperature.ds", overwrite=True)
	hdu.close()
	return hr_lc[:-5] + "temperature.ds"
	