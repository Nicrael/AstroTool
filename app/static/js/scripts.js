document.addEventListener('DOMContentLoaded', function() {
    fetch('/get_settings')
        .then(response => response.json())
        .then(data => {
            const wdBreadcrumbs = data['working_directory'].split('/').pop();
            const dataBreadcrumbs = data['xml_file'].split('/').pop();
            //if no element is present in the DOM dont log an error
            if (document.getElementById('wd_breadcrumbs') === null) {
                return;
            }
            document.getElementById('wd_breadcrumbs').innerText = wdBreadcrumbs;
            if (document.getElementById('data_breadcrumbs') === null) {
                return;
            }
            document.getElementById('data_breadcrumbs').innerText = dataBreadcrumbs;
        })
        .catch(error => console.error('Error fetching settings:', error));
});
