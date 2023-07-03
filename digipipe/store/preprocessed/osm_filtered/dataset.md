# OpenStreetMap gefiltert

OSM data nach bestimmten Tags (s. [config.yml](config.yml) -> `tags`) gefiltert,
zu LAEA Europe (EPSG:3035) umprojiziert und in ein Geopackage konvertiert.

**Achtung:** Konvertierungs- und Extraktionsprozess benötigt 50 GB Speicherplatz
und kann 1 Stunde in Anspruch nehmen.

In einem zweiten Schritt werden diese konvertierten OSM-Daten erneut nach
bestimmten Tags gefiltert und in separate Dateien geschrieben:
- Gebäude: `osm_buildings.gpkg` (s. [config.yml](config.yml) -> `osm_buildings`)
