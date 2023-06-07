# Wärmebedarf

Wärmebedarfe (Endenergie) und -zeitreihen für Haushalte, GHD und Industrie je
Gemeinde.

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor erfolgt anhand
landesweiter Prognosen aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md).

## Haushalte

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme, desaggregiert auf Gemeinden
  mittels Wärmebedarfs-Rasterdaten aus 2015 (Wärmebedarfsdichte 1ha) aus
  [Peta5](../../raw/seenergies_peta5/dataset.md)
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion der
  Gebäudewärmebedarfe aus
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Gasbedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md)-Daten von 2022 die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

## GHD

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme, desaggregiert auf Gemeinden
  mittels Wärmebedarfs-Rasterdaten aus 2015 (Wärmebedarfsdichte 1ha) aus
  [Peta5](../../raw/seenergies_peta5/dataset.md)
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion der
  Gebäudewärmebedarfe aus
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Gasbedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md)-Daten von 2022 die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

## Industrie

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme. Die Desaggregation auf
  Landkreisebene erfolgt anhand des Gesamtenergiebedarfs im verarbeitenden
  Gewerbe aus [Regionalstatistik](../../preprocessed/regiostat/dataset.md).
  Die anschließende Desaggregation auf Gemeindeebene wird mittels
  Beschäftigtenzahlen im verarbeitenden Gewerbe in 2022 aus
  [Regionalstatistik](../../preprocessed/regiostat/dataset.md) vorgenommen.
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  industriellen Gesamtenergiebedarfs aus
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Im Unterschied zu Haushalten und GHD liegen die Daten für den
  Wärme- und Stromanteil nicht getrennt vor, sodass auf den
  Gesamtenergiebedarf zurückgegriffen wird.
  Es wird das Szenario "TN-Strom" als Grundlage für den Status quo verwendet und
  Werte für 2022 interpoliert. Die Zielwerte werden dem Szenario "T45-Strom"
  entnommen.
- Gemittelte, normierte Gasbedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../preprocessed/demandregio/dataset.md)-Daten von 2022 die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.
