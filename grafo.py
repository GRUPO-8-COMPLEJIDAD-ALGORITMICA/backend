import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely.geometry import Point
from math import radians, cos, sin, asin, sqrt
import folium
from sklearn.neighbors import BallTree
import numpy as np

def haversine(lat1, lon1, lat2, lon2):
    R = 6371  # Radio de la tierra en km
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * asin(sqrt(a))
    return R * c

# --- Cargar shapefiles ---
try:
    riesgo = gpd.read_file("puntos_riesgo/points.shp").to_crs(epsg=4326)
    recursos = gpd.read_file("puntos_respuesta/points.shp").to_crs(epsg=4326)
    caminos = gpd.read_file("caminos_lima/roads.shp").to_crs(epsg=4326)
except FileNotFoundError as e:
    print(f"Error al cargar archivos shapefile: {e}")
    exit()

# --- Preparar datos ---
riesgo["tipo"] = "riesgo"
recursos["tipo"] = "recurso"
riesgo = riesgo.loc[:, ~riesgo.columns.duplicated()].reset_index(drop=True)
recursos = recursos.loc[:, ~recursos.columns.duplicated()].reset_index(drop=True)

# Seleccionar un nodo de riesgo específico (por ejemplo, el primero)
nodo_riesgo_index = 0
if nodo_riesgo_index >= len(riesgo):
    print(f"Índice de nodo de riesgo fuera de rango: {nodo_riesgo_index}")
    exit()

nodo_riesgo = riesgo.iloc[nodo_riesgo_index]
riesgo_point = nodo_riesgo.geometry

# --- Buscar recursos cercanos con BallTree ---
recursos_coords = np.array([[geom.y, geom.x] for geom in recursos.geometry])
tree = BallTree(np.radians(recursos_coords), metric='haversine')
radio_km = 5
indices_cercanos = tree.query_radius(np.radians([[riesgo_point.y, riesgo_point.x]]), r=radio_km / 6371.0)[0]

# --- Crear grafo ---
G = nx.Graph()
G.add_node("riesgo", geometry=riesgo_point, tipo="riesgo", nombre=nodo_riesgo.get("nombre", "riesgo"))

# --- Índice espacial de caminos para optimizar ---
caminos_sindex = caminos.sindex

for idx in indices_cercanos:
    recurso = recursos.iloc[idx]
    geom = recurso.geometry
    node_id = f"recurso_{idx}"

    # Buscar solo caminos cercanos al recurso (usando bounding box)
    buffer_geom = geom.buffer(0.005)  # Aproximadamente 500 m
    posibles_idx = list(caminos_sindex.intersection(buffer_geom.bounds))

    conectado = False
    for cidx in posibles_idx:
        camino = caminos.iloc[cidx]
        if geom.distance(camino.geometry) < 0.005:
            dist = haversine(riesgo_point.y, riesgo_point.x, geom.y, geom.x)
            G.add_node(node_id, geometry=geom, tipo="recurso", nombre=recurso.get("nombre", node_id))
            G.add_edge("riesgo", node_id, weight=dist)
            conectado = True
            break  # Solo una conexión válida es suficiente

# --- Crear mapa ---
centro = riesgo_point
m = folium.Map(location=[centro.y, centro.x], zoom_start=13)

# Agregar nodo de riesgo (rojo)
folium.Marker(
    [riesgo_point.y, riesgo_point.x],
    popup=f"Riesgo: {nodo_riesgo.get('nombre', 'riesgo')}",
    icon=folium.Icon(color="red")
).add_to(m)

# Agregar nodos de respuesta (verde si conectados, negro si no)
for i, row in recursos.iterrows():
    color = "green" if f"recurso_{i}" in G.nodes else "black"
    folium.CircleMarker(
        location=[row.geometry.y, row.geometry.x],
        radius=5,
        color=color,
        fill=True,
        fill_color=color,
        fill_opacity=0.7,
        popup=f"Recurso: {row.get('nombre', i)}"
    ).add_to(m)

# Dibujar aristas (azul)
for u, v, data in G.edges(data=True):
    p1 = G.nodes[u]["geometry"]
    p2 = G.nodes[v]["geometry"]
    folium.PolyLine([(p1.y, p1.x), (p2.y, p2.x)], color="blue", weight=2).add_to(m)

# Guardar
m.save("grafo_riesgo_con_recursos.html")
print("Mapa creado: grafo_riesgo_con_recursos.html")
