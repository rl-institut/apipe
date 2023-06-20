# Potenzialgebiete Windenergie

Potenzialgebiete für die Errichtung von Windenergieanlagen, basierend auf den
Teilplänen Wind der Regionalen Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg
aus
[rpg_abw_regional_plan](../../preprocessed/rpg_abw_regional_plan/dataset.md).

Dateien:
- STP Wind 2018 - Vorrang-/Eignungsgebiete:
  `potentialarea_wind_stp_2018_vreg.gpkg`
- STP Wind 2027 - Planabsicht Vorranggebiete:
  `potentialarea_wind_stp_2027_vr.gpkg`
- STP Wind 2027 - Planabsicht Repoweringgebiete:
  `potentialarea_wind_stp_2027_repowering.gpkg`
- STP Wind 2027 - Suchraum Wald:
  `potentialarea_wind_stp_2027_search_area_forest_area.gpkg`
- STP Wind 2027 - Suchraum Offenland:
  `potentialarea_wind_stp_2027_search_area_open_area.gpkg`

Die darin verwendeten Attributtexte werden in die Datei
`potentialarea_wind_attribute_captions.json` exportiert.

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
ha) berechnet, siehe `potentialarea_wind_area_stats_muns.json`. Die Gemeinden
werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.
