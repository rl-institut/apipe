# Technologiedaten

Allgemeine Technologiedaten, Datei: `technology_data.json`

## Jahresvolllaststunden (`full_load_hours`)

Anhand typischer heutiger und prognostizierter Werte für Sachsen-Anhalt werden
folgende Jahresvolllaststunden angenommen:

| Technologie     | Jahr | Volllaststunden | Quelle(n)                                                                                                                                                                                                                                                                                 |
|-----------------|------|----------------:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Windenergie     | 2022 |            1800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wind/auswahl/811-durchschnittliche_ja/#goto_811)                                                                                                                                              |
|                 | 2045 |            2300 | [PV- und Windflächenrechner](https://zenodo.org/record/6794558)                                                                                                                                                                                                                           |
| Freiflächen-PV  | 2022 |             980 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/solar/auswahl/813-durchschnittliche_ja/#goto_813), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/aktuelle-fakten-zur-photovoltaik-in-deutschland.pdf) |
|                 | 2045 |             980 | [PV- und Windflächenrechner](https://zenodo.org/record/6794558), [Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)                                                                                               |
| Aufdach-PV      | 2022 |             910 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/solar/auswahl/813-durchschnittliche_ja/#goto_813), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/aktuelle-fakten-zur-photovoltaik-in-deutschland.pdf) |
|                 | 2045 |             910 | [Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)                                                                                                                                                                                                                                                                                   |
| Laufwasserkraft | 2022 |            3800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wasser/auswahl/840-durchschnittliche_ja/#goto_840)                                                                                                                                            |
|                 | 2045 |            3800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wasser/auswahl/840-durchschnittliche_ja/#goto_840)                                                                                                                                            |

TBD: Generalisieren - automatische Generierung anhand von Global Wind Atlas /
Global Solar Atlas.

## Leistungsdichte (`power_density`)

Installierbare Leistung pro Fläche / spezifischer Flächenbedarf:
- Windenergie: 21 MW/km²
- PV-Freiflächenanlagen: 100 MW/km²
- PV-Aufdachanlagen: 140 MW/km²
- Solarthermie: ? MW/km²

Quelle: [PV- und Windflächenrechner](https://zenodo.org/record/6794558).

## Kosten, Emissionen und Wirkungsgrade

Siehe Datensatz [costs_efficiencies](../costs_efficiencies/dataset.md).
