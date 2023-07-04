# OpenStreetMap Gebäude

OSM Gebäude aus [osm_filtered](../../preprocessed/osm_filtered/dataset.md)
mittels OGR extrahieren und nach Tags (s. [config.yml](config.yml)) filtern.

Ziel ist die Ermittlung des regionalen Anteils Gebäudegrundflächen an der
gesamten Gebäudegrundfläche in Deutschland.

Schritte:
- Extraktion aller Gebäude in Deutschland -> `osm_buildings.gpkg`
- Zentroide und Fläche je Gebäude erstellen -> `osm_buildings_centroids.gpkg`
- Mit Region verschneiden -> `osm_buildings_centroids_region.gpkg`
- Flächensumme berechnen -> `osm_buildings_ground_area_region.gpkg`,
  `osm_buildings_ground_area_country.gpkg`
- Regionalen Anteil berechnen -> `osm_buildings_ground_area_share_region.json`

**Achtung:** Konvertierungs- und Extraktionsprozess benötigt ~15 GB
Speicherplatz und kann viel Zeit in Anspruch nehmen.
