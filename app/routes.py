import os
import json
from flask import Blueprint, render_template, request, redirect, url_for, current_app, send_from_directory, jsonify
from config import config_handler
from data import data_handler

main = Blueprint('main', __name__)

@main.route('/')
def index():
    return render_template('index.html')

@main.route('/observations', methods=['GET', 'POST'])
def observations():
    if request.method == 'POST':
        # Gestisci i dati inviati dal form
        pass
    return render_template('observations.html')

@main.route('/settings', methods=['GET', 'POST'])
def settings():
    if request.method == 'GET':
        return render_template('settings.html')
    if request.method == 'POST':
        folder = request.form.get('working_directory')
        file = request.form.get('xml_file')
        config_handler.set_settings_json(folder, file)
        #TODO:if file does not exist, create it and initialize it
        #TODO:if file exists, check if it is valid
        return redirect(url_for('main.settings'))
    
@main.route('/get_settings', methods=['GET'])
def get_settings():
    settings_json = config_handler.get_settings_json()
    return jsonify(settings_json)

@main.route('/targets', methods=['GET', 'POST'])
def targets():
    if request.method == 'GET':
        return render_template('targets.html')
    if request.method == 'POST':
        pass

@main.route('/get_targets')
def get_targets():
    names = data_handler.get_targets_list()
    return jsonify(names)