# Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft aus Datensatz
[rpg_abw_pv_roof_potential](../../raw/rpg_abw_pv_roof_potential/dataset.md).

Die Gebäudezentroide werden mit den Gemeindegrenzen verschnitten und den
Gemeinden zugeordnet, siehe `potentialarea_pv_roof_area_stats_muns.csv`.
Je Gemeinde wird berechnet:
- Anzahl der Dächer/Gebäude: `roof_count`
- Gesamtgrundfläche: `building_area_sqm`
- Anzahl der Gebäude unter Denkmalschutz: `historic_preservation`
- Installierbare Gesamtleistung je Ausrichtung ("south", "north", "east",
  "west", "flat") in MW: `installable_power_<Ausrichtung>`
- Installierbare Gesamtleistung gesamt in MW: `installable_power_total`
- Erzielbarer Energieertrag je Ausrichtung in MWh ("south", "north", "east",
  "west", "flat"): `energy_annual_<Ausrichtung>`
- Erzielbarer Energieertrag gesamt in MWh: `energy_annual_total`

Des Weiteren wird je Gemeinde der reltive Anteil der bereits installierten
Anlagenleistung an der theoretisch installierbaren Leistung (bei
100% Dachnutzung) berechnet und in `pv_roof_deployment_stats_muns.csv`
geschrieben.

Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.
