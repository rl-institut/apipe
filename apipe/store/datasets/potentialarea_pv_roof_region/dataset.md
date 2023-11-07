# Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft aus Datensatz
[rpg_abw_pv_roof_potential](../../raw/rpg_abw_pv_roof_potential/dataset.md).

Die Gebäudezentroide werden mit den Gemeindegrenzen verschnitten und den
Gemeinden zugeordnet.
Ergebnisdaten:

- Alle Gebäude: `potentialarea_pv_roof_area_stats_muns.csv`
- Alle nicht denkmalgeschützten Gebäude:
  `potentialarea_pv_roof_wo_historic_area_stats_muns.csv`

Des Weiteren wird je Gemeinde der relative Anteil der bereits installierten
Anlagenleistung an der theoretisch installierbaren Leistung (bei
100% Dachnutzung) berechnet.
Ergebnisdaten:

- Alle Gebäude: `potentialarea_pv_roof_deployment_stats_muns.csv`
- Alle nicht denkmalgeschützten Gebäude:
  `potentialarea_pv_roof_wo_historic_deployment_stats_muns.csv`

Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

## Ausbauziele

Es werden PV-Ausbauziele für die Region berechnet, indem die Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der Gebäudegrundflächen disaggregiert werden. Hierzu wird der Anteil der
Gebäudegrundflächen in der Region an der bundesweiten Gebäudegrundflächen
berechnet (s. Datensatz [osm_buildings](../osm_buildings/dataset.md)) und die
Ziele linear skaliert. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen-
und Aufdach-PV unterschieden wird, wird ein Verhältnis von 50:50 angenommen,
d.h. bundesweit 214 GW auf Aufdach-PV entfallen.

Der Anteil beträgt 0,62 % und das Leistungsziel damit 1327 MW, s.
`potentialarea_pv_roof_regionalized_targets.json`.
