# Potenzialgebiete Windenergie

Potenzialgebiete für die Errichtung von Windenergieanlagen, basierend auf den
Teilplänen Wind der Regionalen Planungsgemeinschaft Oderland-Spree aus
[rpg_ols_regional_plan](../../preprocessed/rpg_ols_regional_plan/dataset.md).

Dateien:

- STP Wind 2018 - Eignungsgebiete:
  `potentialarea_wind_stp_2018_vreg.gpkg`
- STP Wind 2024 - Planabsicht Vorranggebiete:
  `potentialarea_wind_stp_2027_vr.gpkg`

Die darin verwendeten Attributtexte werden in die Datei
`potentialarea_wind_attribute_captions.json` exportiert.

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet, siehe `potentialarea_wind_area_stats_muns.csv`. Die Gemeinden
werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.
