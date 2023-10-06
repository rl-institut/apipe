# Technologiedaten

## Jahresvolllaststunden

Anhand typischer heutiger und prognostizierter Werte für Sachsen-Anhalt werden
folgende Jahresvolllaststunden angenommen:

| Technologie     | Jahr | Volllaststunden | Quelle(n) für Annahme                                                                                                                                                                                                                                                                                       | Anmerkung                                                      |
|-----------------|------|----------------:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------|
| Windenergie     | 2022 |            1800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wind/auswahl/811-durchschnittliche_ja/#goto_811)                                                                                                                                                                |                                                                |
|                 | 2045 |            2300 | [PV- und Windflächenrechner](https://zenodo.org/record/6794558)                                                                                                                                                                                                                                             |                                                                |
| Freiflächen-PV  | 2022 |             980 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/solar/auswahl/813-durchschnittliche_ja/#goto_813), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/aktuelle-fakten-zur-photovoltaik-in-deutschland.pdf)                   |                                                                |
|                 | 2045 |             980 | [PV- und Windflächenrechner](https://zenodo.org/record/6794558), [Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)                                                                                                                 |                                                                |
| Aufdach-PV      | 2022 |             910 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/solar/auswahl/813-durchschnittliche_ja/#goto_813), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/aktuelle-fakten-zur-photovoltaik-in-deutschland.pdf)                   |                                                                |
|                 | 2045 |             910 | [Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)                                                                                                                                                                                  |                                                                |
| Laufwasserkraft | 2022 |            3800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wasser/auswahl/840-durchschnittliche_ja/#goto_840)                                                                                                                                                              |                                                                |
|                 | 2045 |            3800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wasser/auswahl/840-durchschnittliche_ja/#goto_840)                                                                                                                                                              |                                                                |
| Bioenergie      | 2022 |            6000 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/bioenergie/auswahl/814-durchschnittliche_ja/#goto_814), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/DE2018_ISE_Studie_Stromgestehungskosten_Erneuerbare_Energien.pdf) | Bioenergie-Stromerzeugung (ohne<br/>biogenen Teil des Abfalls) |
|                 |      |                 |                                                                                                                                                                                                                                                                                                             |                                                                |

Datei: `technology_data.json` --> `full_load_hours`

TBD: Generalisieren - automatische Generierung anhand von Global Wind Atlas /
Global Solar Atlas.

## Leistungsdichte

Installierbare Leistung pro Fläche / spezifischer Flächenbedarf:

- Windenergie: 21 MW/km²
- PV-Freiflächenanlagen: 100 MW/km²
- PV-Aufdachanlagen: 140 MW/km²
- Solarthermie: ? MW/km²

Quelle: [PV- und Windflächenrechner](https://zenodo.org/record/6794558)

Datei: `technology_data.json` --> `power_density`

## Nennleistung Windenergieanlage

Als Zukunftsanlage für 2045 wird eine Enercon E126 6500 (6,5 MW) angenommen.
Diese wird für die Berechnung der Anlagenanzahl in den Ergebnissen
verwendet.

Datei: `technology_data.json` --> `nominal_power_per_unit`

## Batterien

- Kleinbatterien/Heimspeicher: Nennkapazität je installierter PV-Peakleistung
  und Speichernennleistung je installierter Speichernennkapazität aus
  [bnetza_mastr](../../digipipe/store/raw/bnetza_mastr/dataset.md) und
  [HTW](https://solar.htw-berlin.de/wp-content/uploads/HTW-Stromspeicher-Inspektion-2023.pdf).
- Großbatterien: Speichernennleistung je installierter Speichernennkapazität
  aus [bnetza_mastr](../../digipipe/store/raw/bnetza_mastr/dataset.md).

Datei: `technology_data.json` --> `batteries`

## Warmwasserspeicher

- Kleinwärmespeicher (dezentral): Speichernennleistung je installierter
  Speichernennkapazität aus
  [DEA](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage)
- Großwärmespeicher (Fernwärme): Speichernennleistung je installierter
  Speichernennkapazität aus
  [DEA](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage)

Datei: `technology_data.json` --> `hot_water_storages`

## Kosten und Wirkungsgrade

Datei: `raw_costs_efficiencies.csv`
