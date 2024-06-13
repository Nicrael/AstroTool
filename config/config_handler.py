import json

CONFIG_FILE = 'config/config.json'

def load_config(file_path):
    with open(file_path, 'r') as file:
        config = json.load(file)
    return config

def save_config(file_path, config):
    with open(file_path, 'w') as file:
        json.dump(config, file, indent=4)
    return True

def get_settings_json():
    config = load_config(CONFIG_FILE)
    return {
        'working_directory': config.get('working_directory', ''),
        'xml_file': config.get('xml_file', '')
    }

def set_settings_json(folder, file):
    config = load_config(CONFIG_FILE)
    config['working_directory'] = folder
    config['xml_file'] = file
    save_config(CONFIG_FILE, config)