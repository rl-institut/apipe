# Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft aus Datensatz
[rpg_abw_pv_roof_potential](../../raw/rpg_abw_pv_roof_potential/dataset.md).

Die Gebäudezentroide werden mit den Gemeindegrenzen verschnitten und den
Gemeinden zugeordnet, siehe `potentialarea_pv_roof_area_stats_muns.csv`.

Des Weiteren wird je Gemeinde der reltive Anteil der bereits installierten
Anlagenleistung an der theoretisch installierbaren Leistung (bei
100% Dachnutzung) berechnet und in
`potentialarea_pv_roof_deployment_stats_muns.csv` geschrieben.

Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.
