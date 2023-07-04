# Speicheranlagen

Speicheranlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Es wird weiterhin geprüft, ob dem Speicher eine oder mehrere PV-Aufdachanlagen
zugeordnet sind, es wird die Anzahl und Summe der Nettonennleistung berechnet.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Weiterhin erfolgt eine Auswertung der installierten Gesamtleistung je Gemeinde:
- Alle Speicher: `bnetza_mastr_storage_stats_muns.csv`
- Großspeicher (>=100 kWh): `bnetza_mastr_storage_large_stats_muns.csv`
- Kleinspeicher (<100 kWh): `bnetza_mastr_storage_small_stats_muns.csv`

`bnetza_mastr_storage_pv_roof.json` enthält die spezifische Speicherkapazität
sowie spezifische Nennleistung der Speicher (bezogen auf die installierte
Leistung von PV-Aufdachanlagen), aggregiert für gesamte Region, für folgende
Randbedingungen:
- Alle PV-Anlagen: `all_storages`
- PV-Anlagen mit 2..20 kWp sowie Batteriespeicher <20 kWh und <20 kW (kann in
  [config.yml](config.yml) unter `home_storages` konfiguriert werden):
  `home_storages`
