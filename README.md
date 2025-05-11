# Montar la imagen si no existe docker build -t python-name-container .
- docker build -t python-grafo-container .
# Para ejecutar el código introduce este comando en esta carpeta 
- docker run --rm -it -v %cd%:/app python-grafo-container python3 grafo2.py
