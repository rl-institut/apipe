# Strombedarf

Strombedarf für Haushalte, GHD und Industrie auf Gemeindeebene.

Datengrundlage
- Strombedarf 2022: [DemandRegio](../../preprocessed/demandregio/dataset.md)
- Strombedarfsprognosen 2035 und 2045:
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor Haushalte, GHD
und Industrie erfolgt anhand landesweiter Prognosen. Dafür wird der anteilige
Energiebedarf der Region in 2022 am Gesamtbedarf berechnet und dieser unter der
Annahme eines gleichbleibenden regionale Anteils anschließend linear skaliert.
Die Ergebnisse liegen auf NUTS 3-Ebene vor und werden anschließend auf Basis
sektorspezifischer Parameter auf Gemeindeebene desaggregiert (s.u.)

## Haushalte

- Jährlicher Strombedarf je Gemeinde in MWh, von Landkreis- auf Gemeindeebene
  disaggregiert anhand von Bevölkerungsprognosen.
- Gemittelte, normierte Bedarfszeitreihe (auf 1 MWh) aus Daten von 2022 die für
  alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

## GHD

- Jährlicher Strombedarf je Gemeinde in MWh, von Landkreis- auf Gemeindeebene
  disaggregiert anhand von sozialversicherungspflichtig Beschäftigten im Jahr
  2022.
- Gemittelte, normierte Bedarfszeitreihe (auf 1 MWh) aus Daten von 2022 die für
  alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden sowohl
  SLP- als auch branchenspezifische Profile. Aufgrund der geringen Differenzen
  zwischen den Landkreisen werden diese gemittelt. Differenzen zwischen
  verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden und
  Feiertagen und werden daher vernachlässigt.

## Industrie

- Jährlicher Strombedarf je Gemeinde in MWh, von Landkreis- auf Gemeindeebene
  disaggregiert anhand der sozialversicherungspflichtig Beschäftigten im Jahr
  2022.
- Gemittelte, normierte Bedarfszeitreihe (auf 1 MWh) aus Daten von 2022 die für
  alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden sowohl
  SLP- als auch branchenspezifische Profile. Aufgrund der geringen Differenzen
  zwischen den Landkreisen werden diese gemittelt. Differenzen zwischen
  verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden und
  Feiertagen und werden daher vernachlässigt.
