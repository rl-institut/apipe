# Wasserkraftanlagen

Wasserkraftanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_hydro_stats_muns.csv`.
