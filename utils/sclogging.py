import logging as log
import os
import datetime

class CustomFormatter(log.Formatter):
    '''
    Class to define custom log formatting with colors for console output and date format for file output.
    '''
    format = "%(asctime)s - %(levelname)s: %(message)s"
    FORMATS = {
        log.DEBUG: "\033[34m" + format + "\033[0m",  # Blue
        log.INFO: "\033[32m" + format + "\033[0m",  # Green
        log.WARNING: "\033[33m" + format + "\033[0m",  # Yellow
        log.ERROR: "\033[31m" + format + "\033[0m",  # Red
        log.CRITICAL: "\033[41m" + format + "\033[0m",  # Red background
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = log.Formatter(log_fmt, datefmt='%y-%m-%Y %H:%M:%S')
        return formatter.format(record)

# Set basic log level to DEBUG
log.basicConfig(level=log.DEBUG)

# Create log directory if not present
log_dir = 'logs'
os.makedirs(log_dir, exist_ok=True)

# Remove the initial console handler
log.getLogger().handlers = []

# Configure console logging with the desired log level
console_handler = log.StreamHandler()
console_handler.setLevel(log.INFO)
console_handler.setFormatter(CustomFormatter())
log.getLogger().addHandler(console_handler)

# Add date to log file names
current_date = datetime.datetime.now().strftime('%Y-%m-%d')
debug_log_file = os.path.join(log_dir, f'debug_{current_date}.log')
info_log_file = os.path.join(log_dir, f'info_{current_date}.log')

# Configure debug file logging
debug_file_handler = log.FileHandler(debug_log_file)
debug_file_handler.setLevel(log.DEBUG)
debug_file_handler.setFormatter(log.Formatter('%(asctime)s - %(levelname)s: %(message)s', datefmt='%y-%m-%Y %H:%M:%S'))
log.getLogger().addHandler(debug_file_handler)

# Configure info file logging
info_file_handler = log.FileHandler(info_log_file)
info_file_handler.setLevel(log.INFO)
info_file_handler.setFormatter(log.Formatter('%(asctime)s - %(levelname)s: %(message)s', datefmt='%y-%m-%Y %H:%M:%S'))
log.getLogger().addHandler(info_file_handler)
