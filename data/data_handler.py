#custom libs
from utils.sclogging import log
from utils.xmmnewton import *
from config import config_handler
from data.data_classes import *

#astronomical libraries
from astroquery.simbad import Simbad
from astropy.coordinates import Angle
from astropy import coordinates as coord
import astropy.units as u

#general libraries
import os
import xml.dom.minidom
from lxml import etree
import xml.etree.ElementTree as ET

def get_targets_list():
    cwd = os.getcwd()
    working_directory = config_handler.get_settings_json().get('working_directory')
    xml_file = config_handler.get_settings_json().get('xml_file')
    os.chdir(working_directory)
    data = Targets().load(xml_file)
    res = []
    for target in data.targets:
        res.append(target.name)
    os.chdir(cwd)
    return res