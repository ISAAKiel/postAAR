# resultierende Rechtecke in polygone verwandeln
# gleichzeitig Ausrichtung bestimmen
# dafür Umwandlung in Länge/Breite notwendig
# optional: get contextual information from point layer (Datierung von Pfosten, Färbung etc.), 
# dafür resultierendes Rechteck mit Buffer von 1 m versehen, Spatial overlay erstellen und Pfosten gruppieren lassen


# lade rechteck, teile es in Einzelpunkte auf, mess die Entfernung zwischen den Punkten, nehme die längste Strecke und bestimme Ausrichtung

#library(sf)
#sf_airports <- st_as_sf(shp_airports) 
#sf_airports_polygons <- st_polygonize(sf_airports)
#shp_airports <- as(sf_airports_polygons, "Spatial")

# reproject data
#sjer_aoi_WGS84 <- spTransform(sjer_aoi,crs(state_boundary_us))

#gzAzimuth(from, to, type = "snyder_sphere")
