document.addEventListener('DOMContentLoaded', function() {
    fetch('/get_targets')
        .then(response => response.json())
        .then(data => {
            const nameListDiv = document.getElementById('name-list');
            const ul = document.createElement('ul');
            data.forEach(name => {
                const li = document.createElement('li');
                li.textContent = name;
                ul.appendChild(li);
            });
            nameListDiv.appendChild(ul);
        })
        .catch(error => console.error('Error fetching names:', error));
});