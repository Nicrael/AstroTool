# Astro Tool

A web-based tool for managing astronomical data and performing analysis on them.

## Installation

Create a virtual environment and install dependencies:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Running the Application

```bash
python run.py
```

## Project Structure

```
AstroTool/
│
├── app/
│   ├── __init__.py
│   ├── routes.py
│   ├── static/
│   │   ├── ...
│   ├── templates/
│   │   ├── ...
│
├── config/
│   ├── __init__.py
│   ├── settings.py
│   ├── config_handler.py
│
├── data/
│   ├── __init__.py
│   ├── data_loader.py
│   ├── data_processor.py
│
├── utils/
│   ├── __init__.py
│   ├── file_utils.py
│   ├── xml_utils.py
│   ├── astronomy_utils.py
│
├── tests/
│   ├── __init__.py
│   ├── test_routes.py
│   ├── test_config.py
│
├── run.py
├── requirements.txt
└── README.md
```
