# Strombedarf

Nettostrombedarfe und -zeitreihen für Haushalte, GHD und Industrie je Gemeinde.

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor erfolgt anhand
landesweiter Prognosen aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md).

## Haushalte

- Jährlicher Strombedarf je Gemeinde in MWh aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md), von Landkreis- auf
  Gemeindeebene disaggregiert anhand von Bevölkerungsprognosen
  ([STALA ST](../../preprocessed/stala_st_pop_prog/dataset.md)).
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  Strombedarfs (ohne Wärmegewinnung) aus
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Strombedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md)-Daten von 2022, die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

## GHD

- Jährlicher Strombedarf je Gemeinde in MWh aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md), von Landkreis- auf
  Gemeindeebene disaggregiert anhand von sozialversicherungspflichtig
  Beschäftigten im Jahr 2022
  ([BA für Arbeit](../../preprocessed/ba_employment/dataset.md)).
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  Strombedarfs (ohne Wärmegewinnung) aus
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Strombedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md)-Daten von 2022, die
  für alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden
  sowohl SLP- als auch branchenspezifische Profile. Aufgrund der geringen
  Differenzen zwischen den Landkreisen werden diese gemittelt. Differenzen
  zwischen verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden
  und Feiertagen und werden daher vernachlässigt.

## Industrie

- Jährlicher Strombedarf je Gemeinde in MWh aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md), von Landkreis- auf
  Gemeindeebene disaggregiert anhand der Beschäftigten im verarbeitenden
  Gewerbe im Jahr 2022
  ([Regionalstatistik](../../preprocessed/regiostat/dataset.md)).
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  industriellen Gesamtenergiebedarfs aus
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Im Unterschied zu Haushalten und GHD liegen die Daten für den
  Wärme- und Stromanteil nicht getrennt vor, sodass auf den
  Gesamtenergiebedarf zurückgegriffen wird.
  Es wird das Szenario "TN-Strom" als Grundlage für den Status quo verwendet und
  Werte für 2022 interpoliert. Die Zielwerte werden dem Szenario "T45-Strom"
  entnommen.
- Gemittelte, normierte Strombedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md)-Daten von 2022, die
  für alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden
  sowohl SLP- als auch branchenspezifische Profile. Aufgrund der geringen
  Differenzen zwischen den Landkreisen werden diese gemittelt. Differenzen
  zwischen verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden
  und Feiertagen und werden daher vernachlässigt.
