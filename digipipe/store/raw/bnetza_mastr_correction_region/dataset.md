# Marktstammdatenregister Datenkorrektur PV

Überprüfung und manuelle Datenkorrektur der Photovoltaikanlagen aus dem
prozessierten Marktstammdatenregister (Datensatz:
[bnetza_mastr](../bnetza_mastr/dataset.md)).

## Plausibiltätsprüfung

Um grobe Fehler herauszufiltern wird überprüft, ob
- Anlage in Betrieb ist (status = "In Betrieb"),
- Anlage Strom produziert,
- Brutto- und Nettokapazität plausibel sind und
- die Kategorisierung, d.h. Zuordnung eine PV-Anlage zu Freifläche oder Dach,
  plausibel ist (manuelle, visuelle Prüfung von geolokalisierten
  PV-Aufdachanlagen anhand von
  [Orthofotos](https://www.geodatenportal.sachsen-anhalt.de/wss/service/ST_LVermGeo_DOP_WMS_OpenData/guest))

## Dateien

- Korrektur Freiflächenanlagen `bnetza_mastr_pv_ground_region_correction.ods`
- Korrektur Aufdachanlagen `bnetza_mastr_pv_roof_region_correction.ods`

mit Spalten
- _mastr_id_: ID aus dem MaStR
- _reason_: Fehler (wrong_type, wrong_position)
- _wrong_attr_: Fehlerhaftes Attribut
- _correction_: Korrigierter Attributwert (None, wenn Korrektur nicht möglich).
  Korrigierte Geometrien liegen in EPSG:3035 vor.
