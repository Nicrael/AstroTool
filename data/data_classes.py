#custom libs
from utils.sclogging import log
from utils.xmmnewton import *

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

SCHEMA_PATH = '/home/andrea/workspace/GUItool/data/data_structure.xsd'

def validate_data_xml(xml_data, xml_schema=None):
    """
    Validate the data_xml against the schema
    """
    xml_doc = etree.parse(xml_data)

    if xml_schema is None:
        try:
            xsd_path = xml_doc.getroot().get('{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation')
        except:
            log.error("No schema found in the XML file, please provide a schema file!")
            return
        print(os.getcwd())
        xmlschema_doc = etree.parse(xsd_path)
        xmlschema = etree.XMLSchema(xmlschema_doc)
    else:
        xmlschema_doc = etree.parse(xml_schema)
        xmlschema = etree.XMLSchema(xmlschema_doc)

    result = xmlschema.validate(xml_doc)

    if result:
        log.debug("XML file is valid.")
    else:
        log.error("XML file is invalid!")
        error = xmlschema.error_log
        log.error(error.last_error)

    return result

class Targets:
    '''
    Class to store a list of targets.
    '''
    def __init__(self, targets=[], dir_path='', file_name=''):
        if not isinstance(targets, list):
            log.error("targets should be a list!")
            return
        self._targets = targets or []
        if not isinstance(dir_path, str):
            log.error("file_path should be a string!")
            return
        self._dir_path = dir_path or ''
        if dir_path != '':
            os.makedirs(self._dir_path, exist_ok=True)
        if not isinstance(file_name, str):
            log.error("file_name should be a string!")
            return
        self._file_name = file_name or ''
    
    @property
    def targets(self):
        return self._targets
    
    @targets.setter
    def targets(self, targets):
        if not isinstance(targets, list):
            log.error("targets should be a list!")
            return
        self._targets = targets

    def add_target(self, target, overwrite=False):
        if not isinstance(target, Target):
            log.error("Target should be an instance of the Target class!")
            return
        for t in self._targets:
            if t.simbad_name == target.simbad_name:
                if overwrite:
                    self._targets.remove(t)
                else:
                    log.warning(f"The target {target.name} already exists in the list of targets!")
                    return
        if target.dir_path == '':
            target_dir_path = self.dir_path + target.name.replace(' ', '_') + '/'
            target.dir_path = target_dir_path
        self._targets.append(target)

    def add_targets(self, targets, overwrite=False):
        if not isinstance(targets, list):
            log.error("targets should be a list!")
            return
        for target in targets:
            self.add_target(target, overwrite)

    def get_target_by_name(self, target_name):
        for target in self._targets:
            if target.name == target_name:
                return target
        log.warning(f"The target {target_name} does not exist in the list of targets!")
        return None
    
    def remove_target_by_name(self, target_name):
        for target in self._targets:
            if target.name == target_name:
                self._targets.remove(target)
                return
        log.warning(f"The target {target_name} does not exist in the list of targets!")

    def remove_targets(self, target_names):
        for target_name in target_names:
            self.remove_target_by_name(target_name)
            
    def remove_all_targets(self):
        self._targets = []

    @property
    def dir_path(self):
        return self._dir_path
    
    @dir_path.setter
    def dir_path(self, dir_path):
        if not isinstance(dir_path, str):
            log.error("dir_path should be a string!")
            return
        self._dir_path = dir_path

    @property
    def file_name(self):
        return self._file_name
    
    @file_name.setter
    def file_name(self, file_name):
        if not isinstance(file_name, str):
            log.error("file_path should be a string!")
            return
        self._file_name = file_name

    def to_xml(self, schema=SCHEMA_PATH):
        targets = ET.Element("targets")
        for target in self.targets:
            targets.append(target.to_xml())
        targets.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
        targets.set("xsi:noNamespaceSchemaLocation", schema)
        targets.set("dir_path", self.dir_path)
        return targets
    
    @classmethod
    def from_xml(cls, targets_elem, dir_path='', file_name=''):
        targets = []
        dir_path = dir_path
        for target_elem in targets_elem.findall("target"):
            targets.append(Target.from_xml(target_elem))
        return cls(targets=targets, dir_path=dir_path, file_name=file_name)
    
    @classmethod
    def load(cls, file_path):
            log.debug(f"Loading the xml file {file_path}")
            cwd = os.getcwd()
            file = file_path.split("/")[-1]
            os.chdir(file_path.replace(file, ""))
            if not os.path.exists(file):
                os.chdir(cwd)
                log.error(f"The file {file} does not exist!")
                return
            if not validate_data_xml(file):
                log.error("The xml file structure does not match the schema!")
                return
            try:
                tree = ET.parse(file)
                root = tree.getroot()
            except:
                os.chdir(cwd)
                log.error("Could not load the xml file!")
                return
            os.chdir(cwd)
            log.info(f'Data loaded from {file_path} file')
            return cls.from_xml(root, file_path.replace(file, "targets/"), file)
         
    def save(self, filename=''):
        if filename:
            self.file_name = filename
        else:
            if self.file_name:
                filename = self.file_name
            else:
                log.error("No filename provided!")
                return
            
        if self.dir_path:
            path = self.dir_path.replace('targets/', '')
            filename = path + filename
        else:
            self.dir_path = os.getcwd() + '/'

        xml_tree = self.to_xml()
        xml_str = xml.dom.minidom.parseString(ET.tostring(xml_tree)).toprettyxml(indent="    ")
        with open(filename, "w") as file:
            file.write(xml_str)
        log.info(f'Data saved in {filename} file')

class Target:
    '''
    Class to store observations of a single target.
    '''
    def __init__(self, name, simbad_name='', ra='', dec='', observations=[], dir_path='', query_simbad=True):
        self._name = name 
        
        if simbad_name:
            self._simbad_name = simbad_name
        else:
            if query_simbad:
                try:
                    self._simbad_name = Simbad.query_object(self._name)["MAIN_ID"][0]
                except:
                    self._simbad_name = name
                    log.warning("Could not find the simbad name for the target!")
            else:
                self._simbad_name = name
        
        if not ra and query_simbad:
            try:
                self._ra = Angle(Simbad.query_object(self._simbad_name)["RA"][0], unit='hour').degree
            except:
                self._ra = ""
                log.warning("Could not find the RA for the target!")
        else:
            self._ra = ra
        
        if not dec and query_simbad:
            try:
                self._dec = Angle(Simbad.query_object(self._simbad_name)["DEC"][0], unit='degree').degree
            except:
                self._dec = ""
                log.warning("Could not find the DEC for the target!")
        else:
            self._dec = dec
        
        if not isinstance(observations, list):
            log.error("observations should be a list!")
            return
        self._observations = observations or []
        
        if not isinstance(dir_path, str):
            log.error("dir_path should be a string!")
            return
        self._dir_path = dir_path or ''

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        self._name = name

    @property
    def simbad_name(self):
        return self._simbad_name
    
    @simbad_name.setter
    def simbad_name(self, simbad_name):
        self._simbad_name = simbad_name

    @property
    def ra(self):
        return self._ra

    @ra.setter    
    def ra(self, ra):
        self._ra = ra

    @property
    def dec(self):
        return self._dec
    
    @dec.setter
    def dec(self, dec):
        self._dec = dec

    @property
    def coordinates(self):
        return coord.SkyCoord(ra=self.ra, dec=self.dec, unit=(u.degree, u.degree), frame='icrs')
    
    @coordinates.setter
    def coordinates(self, coordinates):
        if not isinstance(coordinates, coord.SkyCoord):
            log.error("Coordinates should be an instance of the SkyCoord class!")
            return
        self._ra = coordinates.ra.degree
        self._dec = coordinates.dec.degree

    @property
    def observations(self):
        return self._observations
    
    @observations.setter
    def observations(self, observations):
        self._observations = observations

    def add_observation(self, observation, overwrite=False):
        if not isinstance(observation, Observation):
            log.error("Observation should be an instance of the Observation class!")
            return
        for obs in self._observations:
            if obs.obsid == observation.obsid and obs.mission == observation.mission:
                if overwrite:
                    self._observations.remove(obs)
                else:
                    log.warning(f"The observation {observation.obsid} already exists in the list of observations!")
                    return
        if observation.dir_path == '':
            observation.dir_path = self.dir_path + observation.obsid + '/'
        self._observations.append(observation)

    def add_observations(self, observations, overwrite=False):
        for observation in observations:
            self.add_observation(observation, overwrite)

    def get_observation_by_obsid(self, obsid):
        for observation in self._observations:
            if observation.obsid == obsid:
                return observation
        log.warning(f"The observation {obsid} does not exist in the list of observations!")
        return None
    
    def remove_observation_by_obsid(self, obsid):
        for observation in self._observations:
            if observation.obsid == obsid:
                self._observations.remove(observation)
                return
        log.warning(f"The observation {obsid} does not exist in the list of observations!")
    
    @property
    def dir_path(self):
        return self._dir_path
    
    @dir_path.setter
    def dir_path(self, dir_path):
        if not isinstance(dir_path, str):
            log.error("dir_path should be a string!")
            return
        self._dir_path = dir_path

    def to_xml(self):
        target = ET.Element("target")
        target.set("name", self.name)
        target.set("simbad_name", self.simbad_name)
        target.set("ra", str(float(self.ra)))
        target.set("dec", str(float(self.dec)))
        target.set("dir_path", self.dir_path)
        if self._observations:
            observations_elem = ET.SubElement(target, "observations")
            for observation in self.observations:
                observations_elem.append(observation.to_xml())
        return target
    
    @classmethod
    def from_xml(cls, target_elem):
        name = target_elem.get("name")
        simbad_name = target_elem.get("simbad_name")
        ra = float(target_elem.get("ra"))
        dec = float(target_elem.get("dec"))
        dir_path = target_elem.get("dir_path")
        observations_elem = target_elem.find("observations")
        if observations_elem:
            observations = [Observation.from_xml(observation_elem) for observation_elem in observations_elem.findall("observation")]
        else:
            observations = []
        res = cls(name, simbad_name, ra, dec, observations, dir_path)
        return res
    
class Observation:
    '''
    Class to store the data of a single observation.
    '''
    def __init__(self,mission, dir_path, obsid, target_name, start_date, end_date, bands, specific_data):
        self._mission = mission
        if not isinstance(obsid, str):
            try:
                log.debug("Converting obsid to string!")
                obsid = str(obsid)
            except:
                log.error("obsid should be a string!")
                return
        self._obsid = obsid
        self._target_name = target_name
        self._start_date = start_date.replace(" ", "T")
        self._end_date = end_date.replace(" ", "T")
        self._bands = bands
        self._specific_data = specific_data

        if not isinstance(dir_path, str):
            log.error("dir_path should be a string!")
            return
        self._dir_path = dir_path or ''

    @property
    def mission(self):
        return self._mission
    
    @mission.setter
    def mission(self, mission):
        self._mission = mission

    @property
    def dir_path(self):
        return self._dir_path
    
    @dir_path.setter
    def dir_path(self, dir_path):   
        if not isinstance(dir_path, str):
            log.error("dir_path should be a string!")
            return
        self._dir_path = dir_path

    @property
    def obsid(self):
        return self._obsid
    
    @obsid.setter
    def obsid(self, obsid):
        if not isinstance(obsid, str):
            log.error("obsid should be a string!")
            return
        self._obsid = obsid

    @property
    def target_name(self):
        return self._target_name
    
    @target_name.setter
    def target_name(self, target_name):
        self._target_name = target_name

    @property
    def start_date(self):
        return self._start_date
    
    @start_date.setter
    def start_date(self, start_date):
        self._start_date = start_date

    @property
    def end_date(self):
        return self._end_date
    
    @end_date.setter
    def end_date(self, end_date):
        self._end_date = end_date

    @property
    def bands(self):
        return self._bands
    
    @bands.setter
    def bands(self, bands):
        self._bands = bands

    @property
    def specific_data(self):
        return self._specific_data
    
    @specific_data.setter
    def specific_data(self, specific_data):
        self._specific_data = specific_data

    @property
    def images(self):
        if self._specific_data:
            return self._specific_data.images
        else:
            return []

    def to_xml(self):
        observation = ET.Element("observation")
        ET.SubElement(observation, "obsid").text = self._obsid
        ET.SubElement(observation, "target_name").text = self._target_name
        ET.SubElement(observation, "start_date").text = self._start_date
        ET.SubElement(observation, "end_date").text = self._end_date
        bands_elem = ET.SubElement(observation, "bands")
        for band in self._bands:
            ET.SubElement(bands_elem, "band").text = band
        if self._specific_data:
            specific_data_elem = ET.SubElement(observation, "specific_data")
            specific_data_elem.append(self._specific_data.to_xml())
        observation.set("mission", self._mission)
        observation.set("dir_path", self._dir_path)
        return observation
    
    @classmethod
    def from_xml(cls, observation_elem):
        mission = observation_elem.get("mission")
        dir_path = observation_elem.get("dir_path")
        obsid = observation_elem.find("obsid").text
        target_name = observation_elem.find("target_name").text
        start_date = observation_elem.find("start_date").text
        end_date = observation_elem.find("end_date").text
        bands = [band.text for band in observation_elem.find("bands").findall("band")]
        specific_data_elem = observation_elem.find("specific_data")
        if specific_data_elem is not None:
            #get first element name
            specific_data_name = specific_data_elem[0].tag
            if specific_data_name == "xmm-newton":
                specific_data = XMMNewtonSpecificData.from_xml(specific_data_elem[0])
            else:
                specific_data = None
        else:
            specific_data = None
        return cls(mission, dir_path, obsid, target_name, start_date, end_date, bands, specific_data)

class SpecificData:
    '''
    Class to store the specific data of an observation.
    '''
    def __init__(self):
        pass

    def from_xml(self, specific_data_elem):
        raise NotImplementedError
    
    def to_xml(self):
        raise NotImplementedError

class XMMNewtonSpecificData(SpecificData):
    '''
    Class to store the specific data of an XMM-Newton observation.
    '''
    def __init__(self, odf=False, pps=False, cif='', sumsas='', instruments=[], reduction_folders=[], images=[]):
        if not isinstance(odf, bool):
            if odf == "true":
                odf = True
            else:
                odf = False
        if not isinstance(pps, bool):
            if pps == "true":
                pps = True
            else:
                pps = False
        self._odf = odf
        self._pps = pps
        self._cif = cif
        self._sumsas = sumsas
        self._instruments = instruments or []
        self._reduction_folders = reduction_folders or []
        self._images = images or []

    @property
    def odf(self):
        return self._odf
    
    @odf.setter
    def odf(self, odf):
        if not isinstance(odf, bool):
            log.error("odf should be a boolean!")
            return
        self._odf = odf

    @property
    def pps(self):
        return self._pps
    
    @pps.setter
    def pps(self, pps):
        if not isinstance(pps, bool):
            log.error("pps should be a boolean!")
            return
        self._pps = pps

    @property
    def cif(self):
        return self._cif
    
    @cif.setter
    def cif(self, cif):
        self._cif = cif
    
    @property
    def sumsas(self):
        return self._sumsas
    
    @sumsas.setter
    def sumsas(self, sumsas):
        self._sumsas = sumsas

    @property
    def instruments(self):
        return self._instruments
    
    @instruments.setter
    def instruments(self, instruments):
        self._instruments = instruments

    @property
    def active_instruments(self):
        return [instrument for instrument in self._instruments if instrument.active]
    
    def add_instrument(self, instrument, overwrite=False):
        if not isinstance(instrument, XMMNewtonInstrument):
            log.error("Instrument should be an instance of the XMMNewtonInstrument class!")
            return
        for i in self._instruments:
            if i.name == instrument.name:
                if overwrite:
                    self._instruments.remove(i)
                else:
                    log.warning(f"The instrument {instrument.name} already exists in the list of instruments!")
                    return
        self._instruments.append(instrument)

    def get_instrument_by_name(self, instrument_name):
        for instrument in self._instruments:
            if instrument.name == instrument_name:
                return instrument
        log.warning(f"The instrument {instrument_name} does not exist in the list of instruments!")
        return None
    
    @property
    def reduction_folders(self):
        return self._reduction_folders
    
    @reduction_folders.setter
    def reduction_folders(self, reduction_folders):
        self._reduction_folders = reduction_folders

    def add_reduction_folder(self, reduction_folder, overwrite=False):
        if not isinstance(reduction_folder, XMMNewtonReduction):
            log.error("Reduction should be an instance of the XMMNewtonReduction class!")
            return
        for r in self._reduction_folders:
            if r.name == reduction_folder.name:
                return
        self._reduction_folders.append(reduction_folder)

    def get_reduction_folders_names(self):
        return [reduction.name for reduction in self._reduction_folders]

    @property
    def images(self):
        return self._images
    
    @images.setter
    def images(self, images):
        self._images = images

    def add_image(self, image, overwrite=False):
        if not isinstance(image, XMMNewtonImage):
            log.error("Image should be an instance of the XMMNewtonImage class!")
            return
        if image.file in [i.file for i in self._images]:
            if overwrite:
                self._images.remove(image)
            else:
                log.warning(f"The image {image.file} already exists in the list of images!")
                return
        self._images.append(image)

    def get_image_by_file(self, file):
        for image in self._images:
            if image.file == file:
                return image
        log.warning(f"The image {file} does not exist in the list of images!")
        return None

    def remove_image_by_file(self, file):
        for image in self._images:
            if image.file == file:
                self._images.remove(image)
                return
        log.warning(f"The image {file} does not exist in the list of images!")

    def to_xml(self):
        xmm_specific_data = ET.Element("xmm-newton")
        xmm_specific_data.set("odf", str(self.odf).lower())
        xmm_specific_data.set("pps", str(self.pps).lower())
        if self.cif != '':
            xmm_specific_data.set("cif", self.cif)
        if self.sumsas != '':
            xmm_specific_data.set("sumsas", self.sumsas)
        if self.instruments != [] and self.instruments is not None:
            instruments_elem = ET.SubElement(xmm_specific_data, "instruments")
            for instrument in self.instruments:
                instruments_elem.append(instrument.to_xml())
        if self.reduction_folders != [] and self.reduction_folders is not None:
            reduction_folders_elem = ET.SubElement(xmm_specific_data, "reduction_folders")
            for reduction_folder in self.reduction_folders:
                reduction_folders_elem.append(reduction_folder.to_xml())
        if self.images != [] and self.images is not None:
            images_elem = ET.SubElement(xmm_specific_data, "images")
            for image in self.images:
                images_elem.append(image.to_xml())
        return xmm_specific_data
    
    @classmethod
    def from_xml(cls, specific_data_elem):
        if specific_data_elem.get("odf") == "true":
            odf = True
        else:
            odf = False
        if specific_data_elem.get("pps") == "true":
            pps = True
        else:
            pps = False
        cif = specific_data_elem.get("cif") if specific_data_elem.get("cif") is not None else ""
        sumsas = specific_data_elem.get("sumsas") if specific_data_elem.get("sumsas") is not None else ""
        if specific_data_elem.find("instruments") is not None:
            instruments = [XMMNewtonInstrument.from_xml(instrument_elem) for instrument_elem in specific_data_elem.find("instruments").findall("instrument")]
        else:
            instruments = []
        if specific_data_elem.find("reduction_folders") is not None:
            reduction_folders = [XMMNewtonReduction.from_xml(reduction_elem) for reduction_elem in specific_data_elem.find("reduction_folders").findall("reduction_folder")]
        else:
            reduction_folders = []
        if specific_data_elem.find("images") is not None:
            images = [XMMNewtonImage.from_xml(image_elem) for image_elem in specific_data_elem.find("images").findall("image")]
        else:
            images = []
        return cls(odf, pps, cif, sumsas, instruments, reduction_folders, images)

class Instrument:
    '''
    Class to store the data of an instrument used in an observation.
    '''
    def __init__(self, name, active):
        self._name = name
        if not isinstance(active, bool):
            if active == "true":
                active = True
            elif active == "false":
                active = False
            else:
                log.error("Active should be a boolean!")
                return
        self._active = active

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        self._name = name

    @property
    def active(self):
        return self._active
    
    @active.setter
    def active(self, active):
        if not isinstance(active, bool):
            log.error("Active should be a boolean!")
            return
        self._active = active

    def to_xml(self):
        instrument_elem = ET.Element("instrument")
        instrument_elem.set("name", self.name)
        instrument_elem.set("active", str(self.active))
        return instrument_elem
    
    @classmethod
    def from_xml(cls, instrument_elem):
        name = instrument_elem.get("name")
        active = instrument_elem.get("active") == "true"
        return cls(name, active)

class XMMNewtonInstrument(Instrument):
    '''
    Class to store the data of an XMM-Newton instrument used in an observation.
    '''
    def __init__(self, name, active, exposure=0, modes=[], filters=[]):
        super().__init__(name, active)
        if isinstance(exposure, int):
            self._exposure = exposure
        else:
            try:
                self._exposure = int(exposure)
            except:
                log.error("Exposure should be a int!")
                return
        self._modes = modes or []
        self._filters = filters or []

    @property
    def exposure(self):
        return self._exposure
    
    @exposure.setter
    def exposure(self, exposure):
        if not isinstance(exposure, int):
            try:
                self._exposure = int(exposure)
            except:
                log.error("Exposure should be a number!")
                return
        else:
            self._exposure = exposure

    @property
    def modes(self):
        return self._modes
    
    @modes.setter
    def modes(self, modes):
        self._modes = modes

    def add_mode(self, mode):
        if not isinstance(mode, str):
            log.error("Mode should be a string!")
            return
        self._modes.append(mode)

    @property
    def filters(self):
        return self._filters
    
    @filters.setter
    def filters(self, filters):
        self._filters = filters

    def add_filter(self, filter_):
        if not isinstance(filter_, str):
            log.error("Filter should be a string!")
            return
        self._filters.append(filter_)

    def to_xml(self):
        instrument_elem = ET.Element("instrument")
        instrument_elem.set("name", self.name)
        instrument_elem.set("active", str(self.active).lower())
        if self.exposure == 0:
            instrument_elem.set("exposure", "0")
            pass
        else:
            instrument_elem.set("exposure", str(self.exposure))
        modes_elem = ET.SubElement(instrument_elem, "modes")
        for mode in self.modes:
            ET.SubElement(modes_elem, "mode").text = mode
        if [filter != "nan" for filter in self.filters] != []:
            filters_elem = ET.SubElement(instrument_elem, "filters")
            for filter in self.filters:
                if filter != "nan" and filter != "":
                    ET.SubElement(filters_elem, "filter").text = filter
        return instrument_elem
    
    @classmethod
    def from_xml(cls, instrument_elem):
        name = instrument_elem.get("name")
        active = instrument_elem.get("active") == "true"
        exposure = instrument_elem.get("exposure")
        if instrument_elem.find("modes") is not None:
            modes = [mode_elem.text for mode_elem in instrument_elem.find("modes").findall("mode")]
        else:
            modes = []
        if instrument_elem.find("filters") is not None:
            filters = [filter_elem.text for filter_elem in instrument_elem.find("filters").findall("filter")]
        else:
            filters = []
        return cls(name, active, exposure, modes, filters)

class Reduction:
    '''
    Class to store the data of a reduction folder.
    '''
    def __init__(self, name, dir_path):
        self._name = name
        if not os.path.exists(dir_path):
            log.warning(f"The folder {dir_path} does not exist!")
            self._dir_path = ''
        self._dir_path = dir_path

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, name):
        self._name = name

    @property
    def dir_path(self):
        return self._dir_path
    
    @dir_path.setter
    def dir_path(self, dir_path):
        if not isinstance(dir_path, str):
            log.error("dir_path should be a string!")
            return
        self._dir_path = dir_path

    def to_xml(self):
        reduction_elem = ET.Element("reduction_folder")
        reduction_elem.set("name", self.name)
        reduction_elem.set("dir_path", self.dir_path)
        return reduction_elem
    
    @classmethod
    def from_xml(cls, reduction_elem):
        name = reduction_elem.get("name")
        dir_path = reduction_elem.get("dir_path")
        return cls(name, dir_path)

class XMMNewtonReduction(Reduction):
    '''
    Class to store the data of an XMM-Newton reduction folder.
    '''
    def __init__(self, name, output_path):
        super().__init__(name, output_path)

class Image:
    '''
    Class to store the data of an image.
    '''
    def __init__(self, file, file_path):
        if not isinstance(file, str):
            log.error("File should be a string!")
            return
        self._file = file
        try:
            with open(file_path, 'r') as f:
                pass
        except:
            log.error(f"The file {file_path} does not exist!")
            return
        self._file_path = file_path
    
    @property
    def file(self):
        return self._file
    
    @file.setter
    def file(self, file):
        if not isinstance(file, str):
            log.error("File should be a string!")
            return
        self._file = file

    @property
    def file_path(self):
        return self._file_path
    
    @file_path.setter
    def file_path(self, file_path):
        try:
            with open(file_path, 'r') as f:
                pass
        except:
            log.error(f"The file {file_path} does not exist!")
            return
        self._file_path = file_path

    def to_xml(self):
        image_elem = ET.Element("image")
        image_elem.set("file", self.file)
        image_elem.set("file_path", self.file_path)
        return image_elem
    
    @classmethod
    def from_xml(cls, image_elem):
        file = image_elem.get("file")
        file_path = image_elem.get("file_path")
        return cls(file, file_path)

class XMMNewtonImage(Image):
    '''
    Class to store the data of an XMM-Newton image.
    '''
    def __init__(self, pipeline, sensor, exp_number, file, file_path, mode, filter, duration):
        super().__init__(file, file_path)
        self._pipeline = pipeline
        self._sensor = sensor
        self._exp_number = exp_number
        self._mode = mode
        self._filter = filter
        self._duration = duration

    @property
    def pipeline(self):
        return self._pipeline
    
    @pipeline.setter
    def pipeline(self, pipeline):
        self._pipeline = pipeline

    @property
    def sensor(self):
        return self._sensor
    
    @sensor.setter
    def sensor(self, sensor):
        self._sensor = sensor

    @property
    def exp_number(self):
        return self._exp_number
    
    @exp_number.setter
    def exp_number(self, exp_number):
        self._exp_number = exp_number

    @property
    def mode(self):
        return self._mode

    @mode.setter    
    def mode(self, mode):
        self._mode = mode

    @property
    def filter(self):
        return self._filter

    @filter.setter    
    def filter(self, filter):
        self._filter = filter   

    @property
    def duration(self):
        return self._duration

    @duration.setter    
    def duration(self, duration):
            self._duration = duration

    def to_xml(self):
        image_elem = ET.Element("image")
        image_elem.set("pipeline", self.pipeline)
        image_elem.set("sensor", self.sensor)
        image_elem.set("exp_number", self.exp_number)
        image_elem.set("file", self.file)
        image_elem.set("file_path", self.file_path)
        image_elem.set("mode", self.mode)
        image_elem.set("filter", self.filter)
        image_elem.set("duration", str(self.duration))
        return image_elem
    
    @classmethod
    def from_xml(cls, image_elem):
        pipeline = image_elem.get("pipeline")
        sensor = image_elem.get("sensor")
        exp_number = image_elem.get("exp_number")
        file = image_elem.get("file")
        file_path = image_elem.get("file_path")
        mode = image_elem.get("mode")
        filter = image_elem.get("filter")
        duration = image_elem.get("duration")
        return cls(pipeline, sensor, exp_number, file, file_path, mode, filter, duration)