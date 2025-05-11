import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely.geometry import Point
from math import radians, cos, sin, asin, sqrt
import folium
from sklearn.neighbors import BallTree
import numpy as np

def haversine(lat1, lon1, lat2, lon2):
    R = 6371  # km
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

# --- Preparar puntos ---
riesgo["tipo"] = "riesgo"
recursos["tipo"] = "recurso"
riesgo = riesgo.loc[:, ~riesgo.columns.duplicated()].reset_index(drop=True)
recursos = recursos.loc[:, ~recursos.columns.duplicated()].reset_index(drop=True)

# Seleccionar un nodo de riesgo estático (por ejemplo, el primer nodo de riesgo)
nodo_riesgo_index = 0  # El índice del nodo de riesgo que queremos usar
if nodo_riesgo_index >= len(riesgo):
    print(f"Error: El índice de nodo de riesgo {nodo_riesgo_index} está fuera de rango.")
    exit()
nodo_riesgo = riesgo.iloc[nodo_riesgo_index]
riesgo_point = nodo_riesgo.geometry

# --- Crear BallTree para búsqueda de vecinos cercanos ---
recursos_coords = np.array([[geom.y, geom.x] for geom in recursos.geometry])
tree_recursos = BallTree(np.radians(recursos_coords), metric='haversine')  # en radianes

# Buscar nodos de respuesta cercanos al nodo de riesgo (radio en km)
radio_busqueda = 5  # 5 kilómetros
indices_cercanos = tree_recursos.query_radius(np.radians([[riesgo_point.y, riesgo_point.x]]), r=radio_busqueda / 6371.0)[0]

# --- Crear grafo ---
G = nx.Graph()
G.add_node(f"riesgo_{nodo_riesgo_index}", geometry=riesgo_point, tipo="riesgo", nombre=nodo_riesgo.get("nombre", f"Riesgo_{nodo_riesgo_index}"))

# --- Conectar nodos de respuesta cercanos ---
nodos_respuesta_cercanos = []
for idx in indices_cercanos:
    recurso_cercano = recursos.iloc[idx]
    node_id = f"recurso_cercano_{idx}"
    G.add_node(node_id, geometry=recurso_cercano.geometry, tipo="recurso_cercano", nombre=recurso_cercano.get("nombre", f"Recurso_cercano_{idx}"))
    nodos_respuesta_cercanos.append(node_id)
    dist = haversine(riesgo_point.y, riesgo_point.x, recurso_cercano.geometry.y, recurso_cercano.geometry.x)
    G.add_edge(f"riesgo_{nodo_riesgo_index}", node_id, weight=dist)

# --- Crear mapa ---
centro = recursos.geometry.unary_union.centroid if not recursos.empty else riesgo_point
m = folium.Map(location=[centro.y, centro.x], zoom_start=10)

# Añadir las zonas de riesgo al mapa (rojo)
folium.GeoJson(riesgo, name="Zonas de Riesgo", style_function=lambda x: {"color": "red", "fillColor": "red", "fillOpacity": 0.2}).add_to(m)

# Añadir los nodos de respuesta
for i, row in recursos.iterrows():
    color = "green" if i in indices_cercanos else "black"
    folium.CircleMarker(
        location=[row['geometry'].y, row['geometry'].x],
        radius=5,
        color=color,
        fill=True,
        fill_color=color,
        fill_opacity=0.7,
        popup=f"Recurso {row.get('nombre', i)}"
    ).add_to(m)

# Dibujar las aristas del grafo (azul)
for u, v, data in G.edges(data=True):
    if u in G.nodes and v in G.nodes:
        p1 = G.nodes[u]['geometry']
        p2 = G.nodes[v]['geometry']
        folium.PolyLine([(p1.y, p1.x), (p2.y, p2.x)], color="blue", weight=1).add_to(m)

# Añadir el nodo de riesgo al mapa (marcador rojo)
folium.Marker(
    [riesgo_point.y, riesgo_point.x],
    popup=f"Riesgo {nodo_riesgo.get('nombre', nodo_riesgo_index)}",
    icon=folium.Icon(color="red")
).add_to(m)

# Guardar el mapa como HTML
m.save("nodos_en_mapa_de_riesgo_recursos_5km.html")
print("Mapa creado con nodos de respuesta cercanos (5km) en verde y lejanos en negro. Archivo: nodos_en_mapa_de_riesgo_recursos_5km.html")