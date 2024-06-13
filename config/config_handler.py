import json

CONFIG_FILE = 'config.json'

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
        'data_dir': config.get('data_dir', ''),
        'xml_data': config.get('xml_data', '')
    }

def set_settings_json(folder, file):
    config = load_config(CONFIG_FILE)
    config['data_dir'] = folder
    config['xml_data'] = file
    save_config(CONFIG_FILE, config)