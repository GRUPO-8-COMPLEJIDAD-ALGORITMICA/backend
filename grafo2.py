import networkx as nx
import geopandas as gpd
import pandas as pd
import folium
import time
from shapely.geometry import Point, LineString
from IPython.display import display
from math import radians, sin, cos, sqrt, atan2

# -------------------
# Función Haversine
# -------------------
def haversine(coord1, coord2):
    R = 6371  # Radio de la Tierra en km
    lat1, lon1 = coord1
    lat2, lon2 = coord2
    dlat = radians(lat2 - lat1)
    dlon = radians(lon2 - lon1)
    a = sin(dlat / 2) ** 2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2) ** 2
    return 2 * R * atan2(sqrt(a), sqrt(1 - a))

# -------------------
# 1. Cargar datos reales
# -------------------
riesgo = gpd.read_file("puntos_riesgo/points.shp").to_crs(epsg=4326)
recursos = gpd.read_file("puntos_respuesta/points.shp").to_crs(epsg=4326)
caminos = gpd.read_file("caminos_lima/roads.shp").to_crs(epsg=4326)

# -------------------
# 2. Construir grafo desde caminos
# -------------------
G = nx.Graph()

# Agregar aristas y nodos desde geometría de líneas
for _, row in caminos.iterrows():
    line: LineString = row.geometry
    coords = list(line.coords)
    for i in range(len(coords) - 1):
        p1 = coords[i]
        p2 = coords[i + 1]
        dist = haversine((p1[1], p1[0]), (p2[1], p2[0]))  # lat, lon
        G.add_edge(p1, p2, dist=dist)

# Guardar posiciones para visualización
positions = {node: (node[1], node[0]) for node in G.nodes}  # (lat, lon)

# -------------------
# 3. Asociar nodos más cercanos dentro de un radio de 5 km
# -------------------
def closest_node(point: Point, radius_km=5):
    nearby_nodes = [n for n in G.nodes if haversine((point.y, point.x), (n[1], n[0])) <= radius_km]
    if nearby_nodes:
        return min(nearby_nodes, key=lambda n: haversine((point.y, point.x), (n[1], n[0])))
    return None  # Si no hay nodos dentro del radio

riesgo['closest_node'] = riesgo.geometry.apply(lambda point: closest_node(point, radius_km=5))
recursos['closest_node'] = recursos.geometry.apply(lambda point: closest_node(point, radius_km=5))

# -------------------
# 4. Correr rutas entre puntos de riesgo y recursos usando Dijkstra
# -------------------
results = []
for _, r in riesgo.iterrows():
    for _, s in recursos.iterrows():
        o = r['closest_node']
        d = s['closest_node']
        if o == d or o is None or d is None:
            continue  # Evitar comparar un nodo consigo mismo o si no hay nodos dentro del radio
        
        try:
            # Calcular la distancia base (haversine) entre los puntos
            baseline = haversine((o[1], o[0]), (d[1], d[0]))

            # Dijkstra
            t0 = time.perf_counter()
            dijkstra_len = nx.shortest_path_length(G, o, d, weight='dist')
            t1 = time.perf_counter()
            time_dij = (t1 - t0) * 1000  # Tiempo en milisegundos

            # A* (con heurística Haversine)
            t0 = time.perf_counter()
            astar_len = nx.astar_path_length(
                G, o, d,
                heuristic=lambda u, v: haversine((u[1], u[0]), (v[1], v[0])),
                weight='dist'
            )
            t1 = time.perf_counter()
            time_astar = (t1 - t0) * 1000  # Tiempo en milisegundos

            # Mejora porcentual
            mejora_dij = 100 * (baseline - dijkstra_len) / baseline
            mejora_astar = 100 * (baseline - astar_len) / baseline

            # Guardar resultados
            results.append({
                'origen': o,
                'destino': d,
                'baseline_km': round(baseline, 2),
                'dijkstra_km': round(dijkstra_len, 2),
                'astar_km': round(astar_len, 2),
                'mejora_dijkstra_%': round(mejora_dij, 2),
                'mejora_astar_%': round(mejora_astar, 2),
                'tiempo_dijkstra_ms': round(time_dij, 2),
                'tiempo_astar_ms': round(time_astar, 2)
            })

        except nx.NetworkXNoPath:
            continue  # Si no hay camino

# Crear DataFrame con los resultados
df = pd.DataFrame(results)
display(df)

# -------------------
# 5. Visualizar en mapa
# -------------------
m = folium.Map(location=[riesgo.geometry.y.mean(), riesgo.geometry.x.mean()], zoom_start=10)

# Nodos
for node in G.nodes:
    lat, lon = node[1], node[0]
    folium.CircleMarker(location=[lat, lon], radius=1, color="gray", fill=True, fill_opacity=0.5).add_to(m)

# Caminos
for u, v in G.edges:
    lat1, lon1 = u[1], u[0]
    lat2, lon2 = v[1], v[0]
    folium.PolyLine(locations=[[lat1, lon1], [lat2, lon2]], color="blue", weight=1).add_to(m)

# Rutas utilizadas
for r in results:
    path = nx.shortest_path(G, r['origen'], r['destino'], weight='dist')
    coords = [(n[1], n[0]) for n in path]  # Convertir de (lon, lat) a (lat, lon)
    folium.PolyLine(locations=coords, color="red", weight=3).add_to(m)

# Puntos de riesgo
for _, row in riesgo.iterrows():
    folium.Marker(
        location=[row.geometry.y, row.geometry.x],
        icon=folium.Icon(color="orange", icon="alert"),
        popup="Riesgo"
    ).add_to(m)

# Puntos de respuesta
for _, row in recursos.iterrows():
    folium.Marker(
        location=[row.geometry.y, row.geometry.x],
        icon=folium.Icon(color="green", icon="plus"),
        popup="Recurso"
    ).add_to(m)

m.save("grafo_riesgo_con_recursos.html")
print("Mapa creado: grafo_riesgo_con_recursos.html")
