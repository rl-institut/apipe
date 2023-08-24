# Wärmebedarf

Wärmebedarfe (Endenergie) Fernwärme und dezentrale Wärme sowie Wärmezeitreihen
für Haushalte, GHD und Industrie je Gemeinde.

## Gesamtwärmebedarf

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor erfolgt anhand
landesweiter Prognosen aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md).

### Haushalte

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme, desaggregiert auf Gemeinden
  mittels Wärmebedarfs-Rasterdaten aus 2015 (Wärmebedarfsdichte 1ha) aus
  [Peta5](../../raw/seenergies_peta5/dataset.md).
  Anm.: Die Desaggregation könnte alternativ über Zensus "Gebäude mit Wohnraum
  nach Heizungsart" (31231-02-01-5, s.
  [regiostat](../../raw/regiostat/dataset.md) erfolgen)
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

### GHD

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

### Industrie

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
- Es erfolgt keine Aufteilung des Wärmebedarfs auf unterschiedliche
  Temperaturniveaus.

## Dezentrale Wärme und Fernwärme

Der Gesamtwärmebedarf wird auf dezentrale Heizsysteme und Fernwärme aufgeteilt.
Fernwärmenetze existieren in Dessau-Roßlau, Bitterfeld-Wolfen, Köthen und
Wittenberg.

Da keine Daten zum tatsächlichen Fernwärmebedarf vorliegen, werden Annahmen auf
Basis folgender Quellen getroffen:

- [Zensus 2011: Gebäude nach Heizungsart](https://www.regionalstatistik.de/genesis//online?operation=table&code=31211-04-01-5-B)
- [BMWK Langfristszenarien: Wärmenachfrage in Wärmenetzen (HH&GHD) (2025)](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/54022/62a2667df6f8c176ff129f7ede944837)
- [STALA ST: Wärmebilanz der Industriebetriebe (2021)](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c256237)
- [STALA ST: Energie- und Wasserversorgung](https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Energie/Berichte/6E403_2020-A.pdf)
- [WindNODE](https://windnode-abw.readthedocs.io/en/latest/energy_system_model.html#district-heating)
- [Peta5: D5 1 District Heating Areas (2020)](https://s-eenergies-open-data-euf.hub.arcgis.com/datasets/b62b8ad79f0e4ae38f032ad6aadb91a0_0/)

Annahmen zu Fernwärmeanteilen (Anteil der Endenergie aus Fernwärme an gesamter
Wärme-Endenergie) je Bedarfssektor:

| Fernwärmenetz     | Haushalte |  GHD | Industrie |
|-------------------|----------:|-----:|----------:|
| Dessau-Roßlau     |      0,36 | 0,36 |      0,19 |
| Bitterfeld-Wolfen |      0,11 | 0,11 |      0,21 |
| Köthen            |      0,07 | 0,07 |      0,21 |
| Wittenberg        |      0,15 | 0,15 |      0,01 |

Die Fernwärmeanteile können in der [config.yml](config.yml) im Abschnitt
`district_heating_share` für jeden Sektor separat angepasst werden. Es wird
vereinfachend angenommen, dass der Anteil an Fernwärme für alle
Szenarien/Zieljahre gleich bleibt.

## Beheizungsstruktur

Die Beheizungsstruktur für 2020 und 2045 wird den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
entnommen (Gebäude: Haushalte und GHD Energiebedarf) und für 2022 interpoliert.
Hierbei wird nach Technologien für dezentrale sowie Fernwärme unterschieden.
Für die Biomasse wird der relative Energiebedarf mit Hilfe von Anteilen der
installierten Leistung von spezifischen Biomasse-Konversionsanlagen
[dbfz_biomasss_capacity_rel](../../preprocessed/dbfz_biomasss_capacity_rel/dataset.md)
je Technologie aufgelöst. Der Vereinfachung halber wird angenommen, dass die
relative installierte Leistung der relativen Energiemenge entspricht.

## Ergebnisdaten

- Haushalte: Wärmebedarf gesamt: `demand_hh_heat_demand.csv`
- Haushalte: Wärmebedarf Fernwärme: `demand_hh_heat_demand_cen.csv`
- Haushalte: Wärmebedarf dezentrale Wärme: `demand_hh_heat_demand_dec.csv`
- Haushalte: Zeitreihen: `demand_hh_heat_timeseries.csv`

- GHD: Wärmebedarf gesamt: `demand_cts_heat_demand.csv`
- GHD: Wärmebedarf Fernwärme: `demand_cts_heat_demand_cen.csv`
- GHD: Wärmebedarf dezentrale Wärme: `demand_cts_heat_demand_dec.csv`
- GHD: Zeitreihen: `demand_cts_heat_timeseries.csv`

- Industrie: Wärmebedarf gesamt: `demand_ind_heat_demand.csv`
- Industrie: Wärmebedarf Fernwärme: `demand_ind_heat_demand_cen.csv`
- Industrie: Wärmebedarf dezentrale Wärme: `demand_ind_heat_demand_dec.csv`
- GHD: Zeitreihen: `demand_ind_heat_timeseries.csv`

- Beheizungsstruktur dezentral (informativ): `demand_heat_structure_dec.csv`
- Beheizungsstruktur zentral (informativ): `demand_heat_structure_cen.csv`
- Beheizungsstruktur dezentral für Weiterverwendung im Energiesystem:
  `demand_heat_structure_esys_dec.csv`
- Beheizungsstruktur Fernwärme für Weiterverwendung im Energiesystem:
  `demand_heat_structure_esys_cen.csv`
