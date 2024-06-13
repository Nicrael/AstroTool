document.addEventListener('DOMContentLoaded', function() {
    fetch('/get_settings')
        .then(response => response.json())
        .then(data => {
            const workingdir = data['working_directory'];
            const xmlfile = data['xml_file'];
            document.getElementById('working_directory').value = workingdir;
            document.getElementById('xml_file').value = xmlfile;
        })
        .catch(error => console.error('Error fetching settings:', error));
});
