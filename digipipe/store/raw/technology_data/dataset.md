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

## Kosten und Wirkungsgrade

Datei: `raw_costs_efficiencies.csv`

#### allgemein

Preise werden aus den Technologie Datenblättern der Danish Energy Agency ([1], [2], [3], [4]) entnommen. Abweichungen werden gesondert genannt.

alle Preise werden auf Euro im Jahr 2020 (dis-)kontiert und damit Inflationsbereinigt.

Für Quellen [1]( https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and), [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants), [3](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage), [4] ist das meist die Umrechnung von 2015 zu 2020. Dafür folgende Formel verwendet:
```
P_(2020) = P_(2015)*f_(infl)
f_(infl) = (1+i_(2015))*(1+i_(2016))...*(1+i_(2019))
f_(infl) = 1,005 * 1,005 * 1.015 * 1,018 * 1,014 = 1,0582
```
[8](https://de.statista.com/themen/112/inflation/#topicOverview)

Werte für 2045 werden durch lineare Interpolation ermittelt.

#### biogas_upgrading plant

Quelle: [4](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-renewable-fuels) "82 Biogas, upgrading"

Aufbereitung von Biogas zu Bio-SNG

#### biogas bpchp_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "06 Gas engines, biogas"

Backpressure Combined heat and power (bpchp) modelliert BHKWs

thermal effiency = electrical_effiency / (c_b+c_v)  (laut [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) S. 390)

#### biogas bpchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "06 Gas engines, biogas"

identische Werte zu biogas bpchp_central. Split fürs Energiesystem, aber eingesetzte technologie identisch

#### biogas_plant

Quelle [4](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-renewable-fuels): "81 Biogas Plant, Basic conf."

Stellt Biogas bereit, welches in KWK (biogas bpchp_central, biogas bpchp_decentral) genutzt werden kann

#### boiler_central

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "44 Natural Gas DH Only"

#### boiler_decentral

Quelle [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants): "202 Gas boiler, ex single", "202 Gas boiler, ex apart", "202 Gas boiler, new single",
"202 Gas boiler, new apart"

Es werden für jedes Szenario jeder Wert aus 4 Komponenten zusammengesetzt.

Diese sind die Kombinationen aus:

- Altbau-Neubau
- Einfamilienhaus-Mehrfamilienhaus

Diese Kompnonten werden durch Faktoren gewichtet zusammengefasst.

Für 2020:

- Verhältnis von Altbau-Neubau aus [7](https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/)
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Neubau aus [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb), verbaute Gasheizungen aggregiert
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Altbau wird als 0.7 / 0.3 angenommen

Für 2045:
- Verhältnis von Altbau-Neubau aus [7](https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/)
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Neubau aus [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb), verbaute Gasheizungen in 2020
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Altbau wird als 0.7 / 0.3 angenommen

volle Berechnungen siehe "boiler_small_script.py" im Code Anhang

#### ch4 bpchp_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "06 Gas engines, natural gas"

Backpressure Combined heat and power (bpchp) modelliert BHKWs

thermal effiency = electrical_effiency / (c_b+c_v)  (laut [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) S. 390)

#### ch4 bpchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "06 Gas engines, natural gas"

identische Werte zu ch4 bpchp_central. Split fürs Energiesystem, aber eingesetzte Technologie identisch

#### ch4 extchp_central

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "05 Gas turb. CC, steam extract., Large", [14] S. 20-21

#### ch4 extchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "05 Gas turb. CC, steam extract., Large"

[14] S. 20-21

identisch wie ch4 extchp_central

#### gt

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "04 Gas turb. simple cycle, L"

gas turbine, offener Prozess

#### heatpump_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "40 Comp. hp, airsource 10 MW"

Wärmepumpentechnologie (Luft-Wasser-WP) aus Langfristigkeitsszenarien

#### heatpump_decentral

Quellen: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants) "207 HP air-water,ex single", "207 HP air-water,ex apart", "207 HP air-water,new single", "207 HP air-water,new apart", "207 HP ground-water,ex single",  "207 HP ground-water,ex apart", "207 HP ground-water,new single", "207 HP ground-water,new apart",
[5], [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb)

Es werden für jedes Szenario jeder Wert aus 8 Komponenten zusammengesetzt.
Diese sind die Kombinationen aus:
- Sole-Umwelt
- Einfamilienhaus-Mehrfamilienhaus (fast alle WP in Einfamilienhäsuern!)
- Altbau-Neubau

Es wird das gemittelte Verhätnis Deutschlandweit der letzten 20 Jahre angenommen (BBSR; Bundesamt für Bauwesen und Raumordnung)

Für 2020 wurden Annahmen für das allgemeine Verhältnis zwischen den Möglichkeiten angenommen:

- Sole-Umwelt sind die aggregierten Absatzzahlen aus [5]
- Einfamilienhaus-Mehrfamilienhaus aus [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb)
- Altbau-Neubau aus [7](https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/)

Mit diesen wird für 2045 wurden Annahmen für das allgemeine Verhältnis zwischen den Möglichkeiten angenommen:
- Sole-Umwelt = 0.87/0.13 (Das sind die Absatzzahlen aus 2022 aus der Branchenstudie)
- Einfamilienhaus-Mehrfamilienhaus = 0.7 / 0.3 (Das ist eine freie Annahme, die eine fortschreitende Verbreitung in Mehrfamilienhäusern annimmt)
- Altbau-Neubau = 0.699 / 0.301 (das gemittelte Verhätnis Deutschlandweit der letzten 20 Jahre)

Die Faktoren in 2045 sind daher:
- Altbau_Umwelt_EFH = 0.4256
- Altbau_Umwelt_MFH = 0.1824
- Altbau_Sole_EFH = 0.0636
- Altbau_Sole_MFH = 0.0272
- Neubau_Umwelt_EFH = 0.1833
- Neubau_Umwelt_MFH = 0.0785
- Neubau_Sole_EFH = 0.0273
- Neubau_Sole_MFH = 0.0117

Berechnung siehe "heatpump_small_script.py" im Code Anhang

#### large_scale_battery

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "180 Lithium Ion Battery", Cebulla [9] S. 181

storage_fixom_cost Berechnung aus UMAS/Oemof_B3 übernommen, ohne Quelle dieser Berechnung gefunden zu haben.

storage_fixom_cost = 0,005 * storage_capacity_cost_overnight

Große Differenzen zwischen Windnode und UMAS, UMAS Methodik übernommen

#### pth_central

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "41 Electric Boilers, small", "41 Electric Boilers, large"

Es wurde ein Mittelwert aus den Electric Biolers small und large gebildet, um relevante Größen in ABW abzubilden.

#### pth_decentral

Quellen: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants): "216 Electric heating,new single", "216 Electric heating,new apart"

Annahmen zu Gebäudebestand siehe heatpump_decentral, nur ohne Kombination mit Altbau, da power to heat in Altbauten vernachlässigbar selten (und wenn in anderen Technologien wie Nachtspeicherheizungen) vorkommt.

Berechnungen siehe "pth_decentral_script" im Code Anhang

#### small_scale_battery

Quelle: [15](https://www.zhb-flensburg.de/fileadmin/content/spezial-einrichtungen/zhb/dokumente/dissertationen/fluri/fluri-2019-wirtschaftlichkeit-dez-stromspeicher.pdf), [17] S. 3

- capacity_cost_overnight: [15] S. 41
- effiency, lost_rate, lifetime: [15] S.v91

#### storage heat_central

Quelle [3](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage): "141 Large hot water tank"

- capacity_cost_overnight und fixom_cost ignoriert, da storage_capacity_cost_overnight, storage_fixom_cost einen Wert hat
- storage heat_decentral

Quelle [3](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage): "141 Large hot water tank"

capacity_cost_overnight und fixom_cost ignoriert, da storage_capacity_cost_overnight, storage_fixom_cost einen Wert hat

Große Differenzen zwischen UMAS und Windnode, UMAS Methodik übernommen

#### hydro ror

Quellen: [16]

- fixom_cost: S. 78
- capacity_cost_overnight: S.75
- lifetime: S. 72

#### lignite oven

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)  "206 Wood stove, single, ex tank"

Der Kohleofen ist eine Komponente, die für die Abbildung des Ist-Zusandes relevant ist.
Die Kohleheizung wird durch gesetzliche Regulierung nicht mehr neu verbaut werden können, wodurch die Komponente für die Optimierung nicht relevant ist.
Auch die Datenlage für die Kohleheizung sehr schlecht ist, die Daten werden daher approximiert.

Keine  direkten Werte vorhanden, daher Modellierung anhand der wood stove Werte

efficiency:

  Differenz der Energie zwischen Holz und Kohle liegt im Heizwert des Brennstoffs. Daher wird die Effizienz der wood stove mit Faktor des Verhältnisses der Heizwerte multipliziert.
  Daten für Heizwerte von BMWK [11](https://www.bmwk.de/Redaktion/DE/Artikel/Energie/energiedaten-gesamtausgabe.html) und [12](https://books.google.de/books?id=n0fVYjrHAlwC&pg=PA58#v=onepage&q&f=false)      ergibt einen Faktor von 4/3

fixom_cost:

  Bestehen großteils aus Brennstoffkosten. Änderung zu wood stove  besteht aus Heizwert (gewonnene Energie pro kg) und Preisdiff pro Kilogramm

  Preise aus brikett-rekord.com [13]

lifetime:

  identisch wie wood stove

  marginal-cost: identisch wie wood stove

  Aus den Annahmen folgt, dass die Investkosten ignoriert werden können.

#### pv_ground

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "22 Utility-scale PV", Vergleich [10]

marginal_cost = 0, da in Quellen nicht vorhanden

Kosten aus [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) im Bereich von [10]

#### pv_rooftop

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "22 PV commercial&industrial rooftop", "22 PV  residential", Vergleich [10]

gewichteter Mittelwert zwischen kommerziellen und Wohnhaus PV.

Gewichtung anhand openMaStR Daten aus Pipeline
```
import geopandas as gpd
import os.path

data_folder = os.path.join("/ROAD/TO/DATA")
data = "bnetza_mastr_pv_roof_region.gpkg"

df = gpd.read_file(os.path.join(data_folder, data))

sum = df[["usage_sector", "status"]].groupby("usage_sector").count().sum()
industrial = (df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Industrie"][0] + \
              df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Sonstige"][0] + \
              df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Landwirtschaft"][0] + \
              df[["usage_sector", "status"]].groupby("usage_sector").count().loc[
                  "Gewerbe, Handel und Dienstleistungen"][0]) \
             / sum
residental = (df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Öffentliches Gebäude"][0] +
              df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Haushalt"][0]) / sum
return [industrial, residental]
```

ergibt 25 % industrial und 75% Haushalte.

marginal_cost = 0, da in Quellen nicht vorhanden

Kosten aus [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) im Bereich von [10]

#### thermalcollector_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) "46 Solar District Heating"

#### thermalcollector_decentral

Quelle: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants) "215 Solar heating,ex single", "215 Solar heating,ex apart", "215 Solar heating,new single", "215 Solar heating,new apart"

Annahmen zu Gebäudebestand siehe heatpump_decentral.

Berechnungen siehe "thermalcollector_decentral_script" im Code Anhang

#### wind onshore

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "20 Onshore turbines", Vergleich [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)

EE Kosten durchweg kleiner als in Windnode in 2020

Windnode bezieht sich auf Frauenhofer ISE aus 2018, Vorgängerstudie zu [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)

Frauenhofer (S. 11) [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html) CAPEX-Range höher als DEA [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) in 2020

1400000-2000000 [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html) zu 1190000 [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and) €/MW

keine Aussagen in Frauenhofer [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html) über 2045

wir wählen DEA als Quelle für die Vergleichbarkeit, da Vergleichbarkeit in der Optimierung der Modellierung Vorrang hat

#### wood extchp_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)  "09a Wood Chips, Medium"

[14] S. 20-21

#### wood extchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)  "09a Wood Chips, Medium"

[14] S. 20-21

identisch zu wood extchp_central

#### wood oven

Quelle: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants), "204 Biomass auto,ex single", "204 Biomass auto,new single", "204 Biomass auto,ex apart", "204 Biomass auto,new apart"

Annahmen zu Gebäudebestand siehe heatpump_decentral.

Berechnungen siehe "wood_oven_script" im Code Anhang

#### Quellen

[1] Danish Energy Agency (2016): "Technology Data - Energy Plants for Electricity and District heating generation", Version 13, von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and

[2] Danish Energy Agency (2016): "Technology Data for heating installations", Version 4, von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants

[3] Danish Energy Agency (2018): "Technology Data – Energy storage", Version 7, von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage

[4] Danish Energy Agency (2017): "Technology Data – Renewable fuels", Versoin 9, von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-renewable-fuels

[5] Karl-Heinz Backhaus (Vaillant), Dr. Hendrik Ehrhardt (Stiebel Eltron), Sven Kersten (NIBE), Steffen
Moser (EnBW), Frank Richert (Wolf), Ingo Rieger (Bosch), Egbert Tippelt (Viessmann), André Jacob
(BWP), Johanna Otting (BWP), Björn Schreinermacher (BWP)(2023): "Branchenstudie 2023: Marktentwicklung – Prognose –Handlungsempfehlungen", Bundesverband Wärmepumpe (BWP) e. V.

[6] Statistisches Landesamt Sachsen-Anhalt: "GENESIS-Tabelle: 31121-0006, Statistik der Baufertigstellungen", von https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb, Stand: 11.04.2023

[7] Statista Research Department(2021): "Struktur des Wohnungsbaus nach Neubau und Sanierung in Deutschland in den Jahren 2001 bis 2020", von https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/, Stand: 03.04.2023 12:26:20

[8] Statista: "Daten und Fakten zur Inflation und den Verbraucherpreisen" , von https://de.statista.com/themen/112/inflation/#topicOverview , Stand: 29.03.2023

[9] Cebulla, Felix (2017): "Storage demand in highly renewable energy scenarios for Europe", OPUS - Online Publikationen der Universität Stuttgart, von https://elib.uni-stuttgart.de/handle/11682/9778

[10] Frauenhofer ISE (2019): "Stromgestehungskosten erneuerbare Energien", von https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html

[11] BMWK (2021): "Energiedaten" von https://www.bmwk.de/Redaktion/DE/Artikel/Energie/energiedaten-gesamtausgabe.html

[12] Michael Herrmann, Jürgen Weber: Öfen und Kamine: Raumheizungen fachgerecht planen und bauen. Beuth Verlag, 201, von  https://books.google.de/books?id=n0fVYjrHAlwC&pg=PA58#v=onepage&q&f=false

[13] www.brikett-rekord.com: "Energiekostenvergleich", von https://www.brikett-rekord.com/de/heizwertvergleich-rekord-briketts.html, letzter Abruf 8.5.2023

[14] WindNode: Modell, Methodik, Daten, ABW; von:  https://wolke.rl-institut.de/apps/files/?dir=/Team_TEO/Projekte/351_Digiplan/Know-How&openfile=316636, letzer Abruf 8.8.2023

[15] Fluri, Verena: "Wirtschaftlichkeit von zukunftsfähigen Geschäftsmodellen dezentraler Stromspeicher" von https://www.zhb-flensburg.de/fileadmin/content/spezial-einrichtungen/zhb/dokumente/dissertationen/fluri/fluri-2019-wirtschaftlichkeit-dez-stromspeicher.pdf, letzter Abruf 8.8.2023

[16] Schröder, Andreas; Kunz, Friedrich; Meiss, Jan; Mendelevitch, Roman; Hirschhausen, Christian von: "Current and Prospective Costs of Electricity Generation until 2050" von https://www.diw.de/documents/publikationen/73/diw_01.c.424566.de/diw_datadoc_2013-068.pdf, letzter Abruf 8.8.2023

[17] Prüggler, Wolfgang (2019): "HEIMSPEICHERSYSTEME UND ZENTRALE BATTERIESPEICHER – KRITISCHE FAKTOREN DER WIRTSCHAFTLICHKEIT" von https://ens.dk/sites/ens.dk/files/Analyser/technology_data_catalogue_for_energy_storage.pdf, letzter Abruf 8.8.2023


### Code Anhang

wood_oven_script.py
```
import pandas as pd
import os.path


def linear_interpolate_2045(wert_1, wert_2):
    zeit_1 = 2040
    zeit_2 = 2050
    wert = wert_1 + (((wert_2 - wert_1) / (zeit_2 - zeit_1)) * (2045 - zeit_1))

    return wert


def get_agg_price_2045(dic):
    # Neubau und Sanierungen allg nach BMI f. Deutschland
    neubau = (0.36 + 0.36 + 0.37 + 0.38 + 0.35 + 0.34 + 0.26 + 0.22 + 0.22 + 0.22 + 0.25 + 0.26 + 0.27 + 0.28 + 0.30 + 0.32 + 0.32 + 0.32 + 0.31 + 0.31) / 20
    altbau = 1 - neubau

    # Verhältnisse Einfamilienhaus-Mehrfamilienhaus nach destatis 2020
    single_new = 693 / 763
    multiple_new = (763 - 693) / 763

    # Einfamilinehaus-Mehrfamilienhaus im Altbau Annahme:
    single_faktor = 0.7
    multiple_faktor = 0.3

    single_new_faktor = neubau * single_new
    multiple_new_faktor = neubau * multiple_new
    single_old_faktor = altbau * single_faktor
    multiple_old_faktor = altbau * multiple_faktor

    single_old = single_old_faktor * dic["single_old_price"]
    multiple_old = multiple_old_faktor * dic["multiple_old_price"]
    single_new = single_new_faktor * dic["single_new_price"]
    multiple_new = multiple_new_faktor * dic["multiple_new_price"]

    preis = single_old + multiple_old + single_new + multiple_new

    return preis


# Daten aus DEA:
# einlesen von Daten
data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

#datensheets
single_old = pd.read_excel(data, "204 Biomass auto,ex single", skiprows=4, nrows=33)
multiple_old = pd.read_excel(data, "204 Biomass auto,ex apart", skiprows=4, nrows=33)
single_new = pd.read_excel(data, "204 Biomass auto,new single", skiprows=4, nrows=33)
multiple_new = pd.read_excel(data, "204 Biomass auto,new apart", skiprows=4, nrows=33)

dic_capacity_cost_overnight_2045 = {
    "single_old_price": linear_interpolate_2045((single_old.iat[19,5]*1000)/(single_old.iat[0,5]/1000), (single_old.iat[19,6]*1000)/(single_old.iat[0,6]/1000)),
    "multiple_old_price": linear_interpolate_2045((multiple_old.iat[19,5]*1000)/(multiple_old.iat[0,5]/1000), (multiple_old.iat[19,6]*1000)/(multiple_old.iat[0,6]/1000)),
    "single_new_price": linear_interpolate_2045((single_new.iat[19,5]*1000)/(single_new.iat[0,5]/1000), (single_new.iat[19,6]*1000)/(single_new.iat[0,6]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[19,5]*1000)/(multiple_new.iat[0,5]/1000), (multiple_new.iat[19,6]*1000)/(multiple_new.iat[0,6]/1000)),
}

dic_effiency_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[3,5], single_old.iat[3,6]),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[3,5], multiple_old.iat[3,6]) ,
    "single_new_price":  linear_interpolate_2045(single_new.iat[3,5], single_new.iat[3,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[3,5], multiple_new.iat[3,6])
}

dic_fixom_cost_2045 = {
    "single_old_price": linear_interpolate_2045((single_old.iat[24,5])/(single_old.iat[0,5]/1000), (single_old.iat[24,6])/(single_old.iat[0,6]/1000)),
    "multiple_old_price":  linear_interpolate_2045((multiple_old.iat[24,5])/(multiple_old.iat[0,5]/1000), (multiple_old.iat[24,6])/(multiple_old.iat[0,6]/1000)),
    "single_new_price": linear_interpolate_2045((single_new.iat[24,5])/(single_new.iat[0,5]/1000), (single_new.iat[24,5])/(single_new.iat[0,5]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[24,5])/(multiple_new.iat[0,5]/1000), (multiple_new.iat[24,6])/(multiple_new.iat[0,6]/1000)),
}
dic_lifetime_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[5,5], single_old.iat[5,6]),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[5,5], multiple_old.iat[5,6]),
    "single_new_price": linear_interpolate_2045(single_new.iat[5,5], single_new.iat[5,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[5,5], multiple_new.iat[5,6]),
}

dic_marginal_cost_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[23,2] / 1000, single_old.iat[23,2] / 1000),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[23,2] / 1000, multiple_old.iat[23,2] / 1000),
    "single_new_price": linear_interpolate_2045(single_new.iat[23,2] / 1000,single_new.iat[23,2] ),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[23,2] / 1000, multiple_new.iat[23,2] / 1000),
}

dic_2045 = [dic_capacity_cost_overnight_2045,dic_effiency_2045, dic_fixom_cost_2045, dic_lifetime_2045, dic_marginal_cost_2045]
val_2045 = []

# Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2045)
```
thermal_collector_small_script.py
```
import pandas as pd
import os.path

def linear_interpolate_2045(wert_1, wert_2):
    zeit_1 = 2040
    zeit_2 = 2050
    wert = wert_1 + (((wert_2 - wert_1) / (zeit_2 - zeit_1)) * (2045 - zeit_1))

    return wert


def get_agg_price_2045(dic):
    # Neubau und Sanierungen allg nach BMI f. Deutschland
    neubau = (0.36 + 0.36 + 0.37 + 0.38 + 0.35 + 0.34 + 0.26 + 0.22 + 0.22 + 0.22 + 0.25 + 0.26 + 0.27 + 0.28 + 0.30 + 0.32 + 0.32 + 0.32 + 0.31 + 0.31) / 20
    altbau = 1 - neubau

    # Verhältnisse Einfamilienhaus-Mehrfamilienhaus nach destatis 2020
    single_new = 693 / 763
    multiple_new = (763 - 693) / 763

    # Einfamilinehaus-Mehrfamilienhaus im Altbau Annahme:
    single_faktor = 0.7
    multiple_faktor = 0.3

    single_new_faktor = neubau * single_new
    multiple_new_faktor = neubau * multiple_new
    single_old_faktor = altbau * single_faktor
    multiple_old_faktor = altbau * multiple_faktor

    single_old = single_old_faktor * dic["single_old_price"]
    multiple_old = multiple_old_faktor * dic["multiple_old_price"]
    single_new = single_new_faktor * dic["single_new_price"]
    multiple_new = multiple_new_faktor * dic["multiple_new_price"]

    preis = single_old + multiple_old + single_new + multiple_new

    return preis


# Daten aus DEA:
# einlesen von Daten
data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

#datensheets
single_old = pd.read_excel(data, "215 Solar heating,ex single", skiprows=4, nrows=33)
multiple_old = pd.read_excel(data, "215 Solar heating,ex apart", skiprows=4, nrows=33)
single_new = pd.read_excel(data, "215 Solar heating,new single", skiprows=4, nrows=33)
multiple_new = pd.read_excel(data, "215 Solar heating,new apart", skiprows=4, nrows=33)

dic_capacity_cost_overnight_2045 = {
    "single_old_price": linear_interpolate_2045((single_old.iat[19,5]*1000)/(single_old.iat[0,5]/1000), (single_old.iat[19,6]*1000)/(single_old.iat[0,6]/1000)),
    "multiple_old_price": linear_interpolate_2045((multiple_old.iat[19,5]*1000)/(multiple_old.iat[0,5]/1000), (multiple_old.iat[19,6]*1000)/(multiple_old.iat[0,6]/1000)),
    "single_new_price": linear_interpolate_2045((single_new.iat[19,5]*1000)/(single_new.iat[0,5]/1000), (single_new.iat[19,6]*1000)/(single_new.iat[0,6]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[19,5]*1000)/(multiple_new.iat[0,5]/1000), (multiple_new.iat[19,6]*1000)/(multiple_new.iat[0,6]/1000)),
}

dic_fixom_cost_2045 = {
    "single_old_price": linear_interpolate_2045((single_old.iat[24,5])/(single_old.iat[0,5]/1000), (single_old.iat[24,6])/(single_old.iat[0,6]/1000)),
    "multiple_old_price":  linear_interpolate_2045((multiple_old.iat[24,5])/(multiple_old.iat[0,5]/1000), (multiple_old.iat[24,6])/(multiple_old.iat[0,6]/1000)),
    "single_new_price": linear_interpolate_2045((single_new.iat[24,5])/(single_new.iat[0,5]/1000), (single_new.iat[24,5])/(single_new.iat[0,5]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[24,5])/(multiple_new.iat[0,5]/1000), (multiple_new.iat[24,6])/(multiple_new.iat[0,6]/1000)),
}
dic_lifetime_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[5,5], single_old.iat[5,6]),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[5,5], multiple_old.iat[5,6]),
    "single_new_price": linear_interpolate_2045(single_new.iat[5,5], single_new.iat[5,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[5,5], multiple_new.iat[5,6]),
}

dic_2045 = [dic_capacity_cost_overnight_2045, dic_fixom_cost_2045, dic_lifetime_2045]
val_2045 = []

# Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2045)
```
pv_rooftop_script.py:
```
import pandas as pd
import geopandas as gpd
import os.path

#trennt residential and industrial rooftop PV nach Nennleistung
def get_proprtion_residential_industrtial(df):
    sum = df[["usage_sector", "status"]].groupby("usage_sector").count().sum()
    industrial = (df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Industrie"][0] + \
                  df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Sonstige"][0] + \
                  df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Landwirtschaft"][0] + \
                  df[["usage_sector", "status"]].groupby("usage_sector").count().loc[
                      "Gewerbe, Handel und Dienstleistungen"][0]) \
                 / sum
    residental = (df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Öffentliches Gebäude"][0] +
                  df[["usage_sector", "status"]].groupby("usage_sector").count().loc["Haushalt"][0]) / sum
    return [industrial, residental]


def get_qgis_df(GeoDataFrame):
    gdf = gpd.read_file(GeoDataFrame, where="geometry_approximated='0'")
    gdf.where(gdf["status"] == "In Betrieb").to_file("bnetza_mastr_pv_roof_region_filtered.gpkg")

def linear_interpolate_2045(wert_1, wert_2):
    zeit_1 = 2040
    zeit_2 = 2050
    wert = wert_1 + (((wert_2 - wert_1) / (zeit_2 - zeit_1)) * (2045 - zeit_1))

    return wert

def get_agg_price_2045(dic, proportion):
    # getting faktoren
    industrial_factor = proportion[0][0]
    residential_factor = proportion[1][0]

    residential = residential_factor * dic["residential_price"]
    industrial = industrial_factor * dic["industrial_price"]

    preis = residential + industrial

    return preis


data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = ["bnetza_mastr_pv_ground_region.gpkg", "bnetza_mastr_pv_roof_region.gpkg"]

df = gpd.read_file(os.path.join(data_folder, data[1]))

# Daten aus DEA:
# einlesen von Daten
data_folder_sheets = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data_sheets = os.path.join(data_folder_sheets, "technology_data_for_el_and_dh.xlsx")

#datensheets
residential = pd.read_excel(data_sheets, "22 Rooftop PV residential", skiprows=4, nrows=42)
industrial = pd.read_excel(data_sheets, "22 Rooftop PV comm.&industrial", skiprows=4, nrows=42)

proportion = get_proprtion_residential_industrtial(df)

dic_capacity_cost_overnight_2045 = {
    "residential_price": linear_interpolate_2045(residential.iat[10,5], residential.iat[10,6])*1000000,
    "industrial_price": linear_interpolate_2045(industrial.iat[10,5], industrial.iat[10,6])*1000000
}
dic_fixom_cost_2045 = {
    "residential_price": linear_interpolate_2045(residential.iat[18,5], residential.iat[18,6]),
    "industrial_price":  linear_interpolate_2045(industrial.iat[18,5], industrial.iat[18,6]),
}

dic_lifetime_2045 = {
    "residential_price": linear_interpolate_2045(residential.iat[3,5], residential.iat[3,6]),
    "industrial_price": linear_interpolate_2045(industrial.iat[3,5], industrial.iat[3,6]),
}

dic_2045 = [dic_capacity_cost_overnight_2045, dic_fixom_cost_2045, dic_lifetime_2045]
val_2045 = []

# Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic, proportion))

print(dic_capacity_cost_overnight_2045, dic_fixom_cost_2045, dic_lifetime_2045)
print(proportion[0][0])
print(val_2045)
```
Pth_decentral_scirpt.py
```
import pandas as pd
import os.path

def linear_interpolate_2045(wert_1, wert_2):
    zeit_1 = 2040
    zeit_2 = 2050
    wert = wert_1 + (((wert_2 - wert_1) / (zeit_2 - zeit_1)) * (2045 - zeit_1))

    return wert


def get_agg_price_2045(dic):
    # Neubau und Sanierungen allg nach BMI f. Deutschland
    neubau = (0.36 + 0.36 + 0.37 + 0.38 + 0.35 + 0.34 + 0.26 + 0.22 + 0.22 + 0.22 + 0.25 + 0.26 + 0.27 + 0.28 + 0.30 + 0.32 + 0.32 + 0.32 + 0.31 + 0.31) / 20

    # Verhältnisse Einfamilienhaus-Mehrfamilienhaus nach destatis 2020
    single_new = 693 / 763
    multiple_new = (763 - 693) / 763

    single_new_faktor = single_new
    multiple_new_faktor = multiple_new


    single_new = single_new_faktor * dic["single_new_price"]
    multiple_new = multiple_new_faktor * dic["multiple_new_price"]

    preis = single_new + multiple_new

    return preis


# Daten aus DEA:
# einlesen von Daten
data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

#datensheets
single_new = pd.read_excel(data, "216 Electric heating,new single", skiprows=4, nrows=33)
multiple_new = pd.read_excel(data, "216 Electric heating,new apart", skiprows=4, nrows=33)

dic_capacity_cost_overnight_2045 = {
    "single_new_price": linear_interpolate_2045((single_new.iat[19,5]*1000)/(single_new.iat[0,5]/1000), (single_new.iat[19,6]*1000)/(single_new.iat[0,6]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[19,5]*1000)/(multiple_new.iat[0,5]/1000), (multiple_new.iat[19,6]*1000)/(multiple_new.iat[0,6]/1000)),
}
dic_effiency_2045 = {
    "single_new_price":  linear_interpolate_2045(single_new.iat[3,5], single_new.iat[3,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[3,5], multiple_new.iat[3,6])
}
dic_fixom_cost_2045 = {
    "single_new_price": linear_interpolate_2045((single_new.iat[24,5])/(single_new.iat[0,5]/1000), (single_new.iat[24,5])/(single_new.iat[0,5]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[24,5])/(multiple_new.iat[0,5]/1000), (multiple_new.iat[24,6])/(multiple_new.iat[0,6]/1000)),
}
dic_lifetime_2045 = {
    "single_new_price": linear_interpolate_2045(single_new.iat[5,5], single_new.iat[5,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[5,5], multiple_new.iat[5,6]),
}
dic_marginal_cost_2045 = {
    "single_new_price": linear_interpolate_2045(single_new.iat[23,2] / 1000,single_new.iat[23,2] ),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[23,2] / 1000, multiple_new.iat[23,2] / 1000),
}

dic_2045 = [dic_capacity_cost_overnight_2045, dic_effiency_2045 , dic_fixom_cost_2045, dic_lifetime_2045, dic_marginal_cost_2045]
val_2045 = []

print(dic_fixom_cost_2045)
# Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2045)
```
heatpumpt_small_script.py:
```
import pandas as pd
import os.path

def get_faktoren_new(df):

    hp_agg = {"single_erdwaerme": 0, "multiple_erdwaerme": 0,  "single_umweltwaerme": 0, "multiple_umweltwaerme": 0,}

    for row in df.itertuples():
        bereich = row[1].split(",")[2]
        energie = row[1].split(",")[3]
        try:
            count_insg = int(row[1].split(",")[4])
            count_single = int(row[1].split(",")[5])
        except:
            ValueError

        if bereich == "Sachsen-Anhalt":
            if energie == "Geothermie":
                hp_agg["single_erdwaerme"] += count_single
                hp_agg["multiple_erdwaerme"] += (count_insg - count_single)
            elif energie == "Umweltthermie (Luft / Wasser)":
                hp_agg["single_umweltwaerme"] += count_single
                hp_agg["multiple_umweltwaerme"] +=  (count_insg - count_single)
            else:
                continue

        else:
            continue

    hp_agg_sum = sum(hp_agg.values())
    air_single_new = hp_agg["single_umweltwaerme"] / hp_agg_sum
    air_multiple_new = hp_agg["multiple_umweltwaerme"] / hp_agg_sum
    ground_single_new = hp_agg["single_erdwaerme"] / hp_agg_sum
    ground_multiple_new = hp_agg["multiple_erdwaerme"] / hp_agg_sum

    return air_single_new, air_multiple_new, ground_single_new, ground_multiple_new

def linear_interpolate_2045(wert_1, wert_2):
    zeit_1 = 2040
    zeit_2 = 2050
    wert = wert_1 + (((wert_2 - wert_1) / (zeit_2 - zeit_1)) * (2045 - zeit_1))

    return wert


def get_agg_price_2020(dic):
    # nach BWP: Absatz von 2010-2020 -> bildet Bestand mit Kosten von Neubau ab
    wp_neubau_abs = 52500+52500+45000+45000+37500+37500+37500+37500+30000+30000+30000
    wp_altbau_abs = 67500+37500+37500+37500+30000+22500+22500+22500+30000+30000+22500
    wp_gesamt_abs = wp_altbau_abs + wp_neubau_abs
    wp_neubau = wp_neubau_abs / wp_gesamt_abs
    wp_geo = 333333 / (333333+750000)
    wp_umwelt = 750000 / (333333+750000) #Umwelt = Luft und Wasser

    # Verhältnisse WP Alt und Neubau in ST nach destatis
        # Daten einlesen
    data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
    hp = os.path.join(data_folder, "2023_04_11_ST_thermische_Primärenergie_neubau_2010-2020.csv")
    df = pd.read_csv(hp, encoding="ISO8859-1", delimiter=";", skiprows=range(0, 10), nrows=2150)

    faktoren_new = get_faktoren_new(df)

    air_water_single_new_faktor = wp_neubau * faktoren_new[0]
    air_water_multiple_new_faktor = wp_neubau * faktoren_new[1]
    ground_water_single_new_fakotr = wp_neubau * faktoren_new[2]
    ground_water_multiple_new_faktor = wp_neubau * faktoren_new[3]

    # wp altbau:
    altbau_air = wp_umwelt - (air_water_single_new_faktor + air_water_multiple_new_faktor)
    altbau_ground = wp_geo - (ground_water_single_new_fakotr + ground_water_multiple_new_faktor)

    # keine Daten, daher wie neubau angenommen (es gibt keinen Grund zu glauben, dass im Mehrfamilien-Altbau mehr WP verbaut werden)
    single_faktor = faktoren_new[0] + faktoren_new[2] # ca 0.95
    multiple_faktor = faktoren_new[1] + faktoren_new[3] # ca 0.05

    air_water_single_old_faktor = altbau_air * single_faktor
    air_water_multiple_old_faktor = altbau_air * multiple_faktor
    ground_water_single_old_fakotr = altbau_ground * single_faktor
    ground_water_multiple_old_faktor = altbau_ground * multiple_faktor


    air_water_single_old = air_water_single_old_faktor * dic["air_water_single_old_price"]
    air_water_multiple_old = air_water_multiple_old_faktor * dic["air_water_multiple_old_price"]
    ground_water_single_old = ground_water_single_old_fakotr * dic["ground_water_single_old_price"]
    ground_water_multiple_old = ground_water_multiple_old_faktor * dic["ground_water_multiple_old_price"]

    air_water_single_new = air_water_single_new_faktor * dic["air_water_single_new_price"]
    air_water_multiple_new = air_water_multiple_new_faktor * dic["air_water_multiple_new_price"]
    ground_water_single_new = ground_water_single_new_fakotr * dic["ground_water_single_new_price"]
    ground_water_multiple_new = ground_water_multiple_new_faktor * dic["ground_water_multiple_new_price"]

    altbau_kosten = air_water_single_old + air_water_multiple_old + ground_water_single_old + ground_water_multiple_old
    neubau_kosten = air_water_single_new + air_water_multiple_new + ground_water_single_new + ground_water_multiple_new

    preis = altbau_kosten + neubau_kosten

    faktoren = [air_water_single_old_faktor,
    air_water_multiple_old_faktor,
    ground_water_single_old_fakotr,
    ground_water_multiple_old_faktor,
    air_water_single_new_faktor,
    air_water_multiple_new_faktor,
    ground_water_single_new_fakotr,
    ground_water_multiple_new_faktor,
    ]

    return preis, faktoren


def get_agg_price_2045(dic):
    # Neubau und Sanierungen allg nach BMI f. Deutschland
    neubau_allg_prozent = (0.36 + 0.36 + 0.37 + 0.38 + 0.35 + 0.34 + 0.26 + 0.22 + 0.22 + 0.22 + 0.25 + 0.26 + 0.27 + 0.28 + 0.30 + 0.32 + 0.32 + 0.32 + 0.31 + 0.31) / 20
    altbau_allg_prozent = 1 - neubau_allg_prozent

    # Sole/Luft nach Absatz 2022 laut BWP
    ground = 0.13
    air = 0.87

    # Einfamilienhaus/Mehrfamilienhaus
    single = 0.7
    multiple = 0.3


    # Faktoren
    air_water_single_old_faktor = altbau_allg_prozent * air * single
    air_water_multiple_old_faktor = altbau_allg_prozent * air * multiple
    ground_water_single_old_fakotr = altbau_allg_prozent * ground * single
    ground_water_multiple_old_faktor = altbau_allg_prozent * ground * multiple
    air_water_single_new_faktor = neubau_allg_prozent * air * single
    air_water_multiple_new_faktor = neubau_allg_prozent * air * multiple
    ground_water_single_new_fakotr = neubau_allg_prozent * ground * single
    ground_water_multiple_new_faktor = neubau_allg_prozent * ground * multiple

    air_water_single_old = air_water_single_old_faktor * dic["air_water_single_old_price"]
    air_water_multiple_old = air_water_multiple_old_faktor * dic["air_water_multiple_old_price"]
    ground_water_single_old = ground_water_single_old_fakotr * dic["ground_water_single_old_price"]
    ground_water_multiple_old = ground_water_multiple_old_faktor * dic["ground_water_multiple_old_price"]
    air_water_single_new = air_water_single_new_faktor * dic["air_water_single_new_price"]
    air_water_multiple_new = air_water_multiple_new_faktor * dic["air_water_multiple_new_price"]
    ground_water_single_new = ground_water_single_new_fakotr * dic["ground_water_single_new_price"]
    ground_water_multiple_new = ground_water_multiple_new_faktor * dic["ground_water_multiple_new_price"]

    altbau_kosten = air_water_single_old + air_water_multiple_old + ground_water_single_old + ground_water_multiple_old
    neubau_kosten = air_water_single_new + air_water_multiple_new + ground_water_single_new + ground_water_multiple_new

    preis = altbau_kosten + neubau_kosten

    faktoren = [air_water_single_old_faktor,
    air_water_multiple_old_faktor,
    ground_water_single_old_fakotr,
    ground_water_multiple_old_faktor,
    air_water_single_new_faktor,
    air_water_multiple_new_faktor,
    ground_water_single_new_fakotr,
    ground_water_multiple_new_faktor,
    ]

    return preis, faktoren


# Daten aus DEA:
# einlesen von Daten
data_folder = os.path.join("/home/local/RL-INSTITUT/aaron.schilling/Dokumente/Projekte/Digipipe")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

#datensheets
air_water_single_old = pd.read_excel(data, "207 HP air-water,ex single", skiprows=4, nrows=33)
air_water_multiple_old = pd.read_excel(data, "207 HP air-water,ex apart", skiprows=4, nrows=33)
ground_water_single_old = pd.read_excel(data, "207 HP air-water,new single", skiprows=4, nrows=33)
ground_water_multiple_old = pd.read_excel(data,  "207 HP air-water,new apart", skiprows=4, nrows=33)
air_water_single_new = pd.read_excel(data, "207 HP ground-water,ex single", skiprows=4, nrows=33)
air_water_multiple_new = pd.read_excel(data, "207 HP ground-water,ex apart", skiprows=4, nrows=33)
ground_water_single_new = pd.read_excel(data, "207 HP ground-water,new single", skiprows=4, nrows=33)
ground_water_multiple_new = pd.read_excel(data, "207 HP ground-water,new apart", skiprows=4, nrows=33)

dic_capacity_cost_overnight_2020 = {
        "air_water_single_old_price": air_water_single_old.iat[19,2]*1000/(air_water_single_old.iat[0,2]/1000),
        "air_water_multiple_old_price": air_water_multiple_old.iat[19,2]*1000/(air_water_multiple_old.iat[0,2]/1000),
        "ground_water_single_old_price": ground_water_single_old.iat[19,2]*1000/(ground_water_single_old.iat[0,2]/1000),
        "ground_water_multiple_old_price": ground_water_multiple_old.iat[19,2]*1000/(ground_water_multiple_old.iat[0,2]/1000),
        "air_water_single_new_price": air_water_single_new.iat[19,2]*1000/(air_water_single_new.iat[0,2]/1000),
        "air_water_multiple_new_price": air_water_multiple_new.iat[19,2]*1000/(air_water_multiple_new.iat[0,2]/1000),
        "ground_water_single_new_price": ground_water_single_new.iat[19,2]*1000/(ground_water_single_new.iat[0,2]/1000),
        "ground_water_multiple_new_price": ground_water_multiple_new.iat[19,2]*1000/(ground_water_multiple_new.iat[0,2]/1000),
}
dic_fixom_cost_2020 = {
        "air_water_single_old_price": air_water_single_old.iat[24,2]/(air_water_single_old.iat[0,2]/1000),
        "air_water_multiple_old_price": air_water_multiple_old.iat[24,2]/(air_water_multiple_old.iat[0,2]/1000),
        "ground_water_single_old_price": ground_water_single_old.iat[24,2]/(ground_water_single_old.iat[0,2]/1000),
        "ground_water_multiple_old_price": ground_water_multiple_old.iat[24,2]/(ground_water_multiple_old.iat[0,2]/1000),
        "air_water_single_new_price": air_water_single_new.iat[24,2]/(air_water_single_new.iat[0,2]/1000),
        "air_water_multiple_new_price": air_water_multiple_new.iat[24,2]/(air_water_multiple_new.iat[0,2]/1000),
        "ground_water_single_new_price": ground_water_single_new.iat[24,2]/(ground_water_single_new.iat[0,2]/1000),
        "ground_water_multiple_new_price": ground_water_multiple_new.iat[24,2]/(ground_water_multiple_new.iat[0,2]/1000),
}
dic_lifetime_2020 = {
        "air_water_single_old_price": air_water_single_old.iat[5,2],
        "air_water_multiple_old_price": air_water_multiple_old.iat[5,2],
        "ground_water_single_old_price": ground_water_single_old.iat[5,2],
        "ground_water_multiple_old_price": ground_water_multiple_old.iat[5,2],
        "air_water_single_new_price": air_water_single_new.iat[5,2],
        "air_water_multiple_new_price": air_water_multiple_new.iat[5,2],
        "ground_water_single_new_price": ground_water_single_new.iat[5,2],
        "ground_water_multiple_new_price": ground_water_multiple_new.iat[5,2],
}
dic_marginal_cost_2020 = {
        "air_water_single_old_price": air_water_single_old.iat[23,2] / 1000,
        "air_water_multiple_old_price": air_water_multiple_old.iat[23,2] / 1000,
        "ground_water_single_old_price": ground_water_single_old.iat[23,2] / 1000,
        "ground_water_multiple_old_price": ground_water_multiple_old.iat[23,2] / 1000,
        "air_water_single_new_price": air_water_single_new.iat[23,2] / 1000,
        "air_water_multiple_new_price": air_water_multiple_new.iat[23,2] / 1000,
        "ground_water_single_new_price": ground_water_single_new.iat[23,2] / 1000,
        "ground_water_multiple_new_price": ground_water_multiple_new.iat[23,2] / 1000,
}
dic_capacity_cost_overnight_2045 = {
        "air_water_single_old_price": linear_interpolate_2045(air_water_single_old.iat[19,5]*1000/(air_water_single_old.iat[0,5]/1000), air_water_single_old.iat[19,6]*1000/(air_water_single_old.iat[0,6]/1000)),
        "air_water_multiple_old_price": linear_interpolate_2045(air_water_multiple_old.iat[19,5]*1000/(air_water_multiple_old.iat[0,5]/1000), air_water_multiple_old.iat[19,6]*1000/(air_water_multiple_old.iat[0,6]/1000)),
        "ground_water_single_old_price": linear_interpolate_2045(ground_water_single_old.iat[19,5]*1000/(ground_water_single_old.iat[0,5]/1000),ground_water_single_old.iat[19,6]*1000/(ground_water_single_old.iat[0,6]/1000)),
        "ground_water_multiple_old_price": linear_interpolate_2045(ground_water_multiple_old.iat[19,5]*1000/(ground_water_multiple_old.iat[0,5]/1000), ground_water_multiple_old.iat[19,6]*1000/(ground_water_multiple_old.iat[0,6]/1000)),
        "air_water_single_new_price": linear_interpolate_2045(air_water_single_new.iat[19,5]*1000/(air_water_single_new.iat[0,5]/1000), air_water_single_new.iat[19,6]*1000/(air_water_single_new.iat[0,6]/1000)),
        "air_water_multiple_new_price": linear_interpolate_2045(air_water_multiple_new.iat[19,5]*1000/(air_water_multiple_new.iat[0,5]/1000),air_water_multiple_new.iat[19,6]*1000/(air_water_multiple_new.iat[0,6]/1000)),
        "ground_water_single_new_price": linear_interpolate_2045(ground_water_single_new.iat[19,5]*1000/(ground_water_single_new.iat[0,5]/1000), ground_water_single_new.iat[19,6]*1000/(ground_water_single_new.iat[0,6]/1000)),
        "ground_water_multiple_new_price": linear_interpolate_2045(ground_water_multiple_new.iat[19,5]*1000/(ground_water_multiple_new.iat[0,5]/1000), ground_water_multiple_new.iat[19,6]*1000/(ground_water_multiple_new.iat[0,6]/1000)),
}
dic_fixom_cost_2045 = {
        "air_water_single_old_price": linear_interpolate_2045(air_water_single_old.iat[24,5]/(air_water_single_old.iat[0,5]/1000), air_water_single_old.iat[24,6]/(air_water_single_old.iat[0,6]/1000)),
        "air_water_multiple_old_price": linear_interpolate_2045(air_water_multiple_old.iat[24,5]/(air_water_multiple_old.iat[0,5]/1000), air_water_multiple_old.iat[24,6]/(air_water_multiple_old.iat[0,6]/1000)),
        "ground_water_single_old_price": linear_interpolate_2045(ground_water_single_old.iat[24,5]/(ground_water_single_old.iat[0,5]/1000), ground_water_single_old.iat[24,6]/(ground_water_single_old.iat[0,6]/1000)),
        "ground_water_multiple_old_price": linear_interpolate_2045(ground_water_multiple_old.iat[24,5]/(ground_water_multiple_old.iat[0,5]/1000), ground_water_multiple_old.iat[24,6]/(ground_water_multiple_old.iat[0,6]/1000)),
        "air_water_single_new_price": linear_interpolate_2045(air_water_single_new.iat[24,5]/(air_water_single_new.iat[0,5]/1000), air_water_single_new.iat[24,6]/(air_water_single_new.iat[0,6]/1000)),
        "air_water_multiple_new_price": linear_interpolate_2045(air_water_multiple_new.iat[24,5]/(air_water_multiple_new.iat[0,5]/1000), air_water_multiple_new.iat[24,6]/(air_water_multiple_new.iat[0,6]/1000)),
        "ground_water_single_new_price": linear_interpolate_2045(ground_water_single_new.iat[24,5]/(ground_water_single_new.iat[0,5]/1000), ground_water_single_new.iat[24,6]/(ground_water_single_new.iat[0,6]/1000)),
        "ground_water_multiple_new_price": linear_interpolate_2045(ground_water_multiple_new.iat[24,5]/(ground_water_multiple_new.iat[0,5]/1000), ground_water_multiple_new.iat[24,6]/(ground_water_multiple_new.iat[0,6]/1000)),
}
dic_lifetime_2045 = {
        "air_water_single_old_price": linear_interpolate_2045(air_water_single_old.iat[5,5], air_water_single_old.iat[5,6]),
        "air_water_multiple_old_price": linear_interpolate_2045(air_water_multiple_old.iat[5,5], air_water_multiple_old.iat[5,6]),
        "ground_water_single_old_price": linear_interpolate_2045(ground_water_single_old.iat[5,5], ground_water_single_old.iat[5,6]),
        "ground_water_multiple_old_price": linear_interpolate_2045(ground_water_multiple_old.iat[5,5], ground_water_multiple_old.iat[5,6]),
        "air_water_single_new_price": linear_interpolate_2045(air_water_single_new.iat[5,5], air_water_single_new.iat[5,6]),
        "air_water_multiple_new_price": linear_interpolate_2045(air_water_multiple_new.iat[5,5], air_water_multiple_new.iat[5,6]),
        "ground_water_single_new_price": linear_interpolate_2045(ground_water_single_new.iat[5,5], ground_water_single_new.iat[5,6]),
        "ground_water_multiple_new_price": linear_interpolate_2045(ground_water_multiple_new.iat[5,5], ground_water_multiple_new.iat[5,6]),
}
dic_marginal_cost_2045 = {
        "air_water_single_old_price": linear_interpolate_2045(air_water_single_old.iat[23,5] / 1000, air_water_single_old.iat[23,6] / 1000),
        "air_water_multiple_old_price": linear_interpolate_2045(air_water_multiple_old.iat[23,5] / 1000, air_water_multiple_old.iat[23,6] / 1000),
        "ground_water_single_old_price": linear_interpolate_2045(ground_water_single_old.iat[23,5] / 1000, ground_water_single_old.iat[23,6] / 1000),
        "ground_water_multiple_old_price": linear_interpolate_2045(ground_water_multiple_old.iat[23,5] / 1000, ground_water_multiple_old.iat[23,6] / 1000),
        "air_water_single_new_price": linear_interpolate_2045(air_water_single_new.iat[23,5] / 1000, air_water_single_new.iat[23,6] / 1000),
        "air_water_multiple_new_price": linear_interpolate_2045(air_water_multiple_new.iat[23,5] / 1000, air_water_multiple_new.iat[23,6] / 1000),
        "ground_water_single_new_price": linear_interpolate_2045(ground_water_single_new.iat[23,5] / 1000, ground_water_single_new.iat[23,6] / 1000),
        "ground_water_multiple_new_price": linear_interpolate_2045(ground_water_multiple_new.iat[23,5] / 1000, ground_water_multiple_new.iat[23,6] / 1000),
}

dic_2020 = [dic_capacity_cost_overnight_2020, dic_fixom_cost_2020, dic_lifetime_2020, dic_marginal_cost_2020]
dic_2045 = [dic_capacity_cost_overnight_2045, dic_fixom_cost_2045, dic_lifetime_2045, dic_marginal_cost_2045]

val_2020 = []
val_2045 = []
faktoren_2020 = get_agg_price_2020(dic_fixom_cost_2045)[1]
faktoren_2045 = get_agg_price_2045(dic_fixom_cost_2045)[1]

# Berechnungen
for dic in dic_2020:
    val_2020.append(get_agg_price_2020(dic)[0])

for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic)[0])

print(val_2020, val_2045)
print(faktoren_2020, faktoren_2045)
```
boiler_small_scirpt.py:
```
import pandas as pd
import os.path

def get_faktoren_new(df):

    gas_agg = {"single": 0, "multiple": 0}

    for row in df.itertuples():
        bereich = row[1].split(",")[2]
        energie = row[1].split(",")[3]
        try:
            count_insg = int(row[1].split(",")[4])
            count_single = int(row[1].split(",")[5])
        except:
            ValueError

        if bereich == "Sachsen-Anhalt":
            if energie == "Gas":
                gas_agg["single"] += count_single
                gas_agg["multiple"] += (count_insg - count_single)
            else:
                continue

        else:
            continue

    gas_agg_sum = sum(gas_agg.values())
    single_new = gas_agg["single"] / gas_agg_sum
    multiple_new = gas_agg["multiple"] / gas_agg_sum

    return single_new, multiple_new,


def linear_interpolate_2045(wert_1, wert_2):
    zeit_1 = 2040
    zeit_2 = 2050
    wert = wert_1 + (((wert_2 - wert_1) / (zeit_2 - zeit_1)) * (2045 - zeit_1))

    return wert


def get_agg_price_2020(dic):
    # Neubau und Sanierungen allg nach BMI f. Deutschland
    neubau = (0.36 + 0.36 + 0.37 + 0.38 + 0.35 + 0.34 + 0.26 + 0.22 + 0.22 + 0.22 + 0.25 + 0.26 + 0.27 + 0.28 + 0.30 + 0.32 + 0.32 + 0.32 + 0.31 + 0.31) / 20
    altbau = 1 - neubau

    # Verhältnisse Einfamilienhaus-Mehrfamilienhaus nach destatis
        # Daten einlesen
    data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
    hp = os.path.join(data_folder, "2023_04_11_ST_thermische_Primärenergie_neubau_2010-2020.csv")
    df = pd.read_csv(hp, encoding="ISO8859-1", delimiter=";", skiprows=range(0, 10), nrows=2150)

    faktoren_new = get_faktoren_new(df)

    # Einfamilinehaus-Mehrfamilienhaus im Altbau Annahme:
    single_faktor = 0.7
    multiple_faktor = 0.3

    single_new_faktor = neubau * faktoren_new[0]
    multiple_new_faktor = neubau * faktoren_new[1]
    single_old_faktor = altbau * single_faktor
    multiple_old_faktor = altbau * multiple_faktor

    single_old = single_old_faktor * dic["single_old_price"]
    multiple_old = multiple_old_faktor * dic["multiple_old_price"]
    single_new = single_new_faktor * dic["single_new_price"]
    multiple_new = multiple_new_faktor * dic["multiple_new_price"]

    preis = single_old + multiple_old + single_new + multiple_new


    return preis


def get_agg_price_2045(dic):
    # Neubau und Sanierungen allg nach BMI f. Deutschland
    neubau = (0.36 + 0.36 + 0.37 + 0.38 + 0.35 + 0.34 + 0.26 + 0.22 + 0.22 + 0.22 + 0.25 + 0.26 + 0.27 + 0.28 + 0.30 + 0.32 + 0.32 + 0.32 + 0.31 + 0.31) / 20
    altbau = 1 - neubau

    # Verhältnisse Einfamilienhaus-Mehrfamilienhaus nach destatis 2020
    gas_single_new = 693 / 763
    gas_multiple_new = (763 - 693) / 763

    # Einfamilinehaus-Mehrfamilienhaus im Altbau Annahme:
    single_faktor = 0.7
    multiple_faktor = 0.3

    single_new_faktor = neubau * gas_single_new
    multiple_new_faktor = neubau * gas_multiple_new
    single_old_faktor = altbau * single_faktor
    multiple_old_faktor = altbau * multiple_faktor

    single_old = single_old_faktor * dic["single_old_price"]
    multiple_old = multiple_old_faktor * dic["multiple_old_price"]
    single_new = single_new_faktor * dic["single_new_price"]
    multiple_new = multiple_new_faktor * dic["multiple_new_price"]

    preis = single_old + multiple_old + single_new + multiple_new

    return preis


# Daten aus DEA:
# einlesen von Daten
data_folder = os.path.join("/ROAD/TO/DATA")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

#datensheets
single_old = pd.read_excel(data, "202 Gas boiler, ex single", skiprows=4, nrows=33)
multiple_old = pd.read_excel(data, "202 Gas boiler, ex apart", skiprows=4, nrows=33)
single_new = pd.read_excel(data, "202 Gas boiler, new single", skiprows=4, nrows=33)
multiple_new = pd.read_excel(data, "202 Gas boiler, new apart", skiprows=4, nrows=33)


dic_capacity_cost_overnight_2020 = {
    "single_old_price": (single_old.iat[19,2]*1000)/(single_old.iat[0,2]/1000),
    "multiple_old_price": (multiple_old.iat[19,2]*100)/(multiple_old.iat[0,2]/1000),
    "single_new_price": (single_new.iat[19,2]*1000)/(single_new.iat[0,2]/1000),
    "multiple_new_price": (multiple_new.iat[19,2]*1000)/(multiple_new.iat[0,2]/1000),
}
dic_effiency_2020 = {
    "single_old_price": single_old.iat[3,2],
    "multiple_old_price": multiple_old.iat[3,2],
    "single_new_price": single_new.iat[3,2],
    "multiple_new_price": multiple_new.iat[3,2],
}
dic_fixom_cost_2020 = {
    "single_old_price": single_old.iat[24,2]/(single_old.iat[0,2]/1000),
    "multiple_old_price": multiple_old.iat[24,2]/(multiple_old.iat[0,2]/1000),
    "single_new_price": single_new.iat[24,2]/(single_new.iat[0,2]/1000),
    "multiple_new_price": multiple_new.iat[24,2]/(multiple_new.iat[0,2]/1000),
}
dic_lifetime_2020 = {
    "single_old_price": single_old.iat[5,2],
    "multiple_old_price": multiple_old.iat[5,2],
    "single_new_price": single_new.iat[5,2],
    "multiple_new_price": multiple_new.iat[5,2],
}
dic_marginal_cost_2020 = {
    "single_old_price": single_old.iat[23,2] / 1000,
    "multiple_old_price": multiple_old.iat[23,2] / 1000,
    "single_new_price": single_new.iat[23,2] / 1000,
    "multiple_new_price": multiple_new.iat[23,2] / 1000,
}

dic_capacity_cost_overnight_2045 = {
    "single_old_price": linear_interpolate_2045((single_old.iat[19,5]*1000)/(single_old.iat[0,5]/1000), (single_old.iat[19,6]*1000)/(single_old.iat[0,6]/1000)),
    "multiple_old_price": linear_interpolate_2045((multiple_old.iat[19,5]*1000)/(multiple_old.iat[0,5]/1000), (multiple_old.iat[19,6]*1000)/(multiple_old.iat[0,6]/1000)),
    "single_new_price": linear_interpolate_2045((single_new.iat[19,5]*1000)/(single_new.iat[0,5]/1000), (single_new.iat[19,6]*1000)/(single_new.iat[0,6]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[19,5]*1000)/(multiple_new.iat[0,5]/1000), (multiple_new.iat[19,6]*1000)/(multiple_new.iat[0,6]/1000)),
}

dic_effiency_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[3,5], single_old.iat[3,6]),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[3,5], multiple_old.iat[3,6]) ,
    "single_new_price":  linear_interpolate_2045(single_new.iat[3,5], single_new.iat[3,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[3,5], multiple_new.iat[3,6])
}

dic_fixom_cost_2045 = {
    "single_old_price": linear_interpolate_2045((single_old.iat[24,5])/(single_old.iat[0,5]/1000), (single_old.iat[24,6])/(single_old.iat[0,6]/1000)),
    "multiple_old_price":  linear_interpolate_2045((multiple_old.iat[24,5])/(multiple_old.iat[0,5]/1000), (multiple_old.iat[24,6])/(multiple_old.iat[0,6]/1000)),
    "single_new_price": linear_interpolate_2045((single_new.iat[24,5])/(single_new.iat[0,5]/1000), (single_new.iat[24,5])/(single_new.iat[0,5]/1000)),
    "multiple_new_price": linear_interpolate_2045((multiple_new.iat[24,5])/(multiple_new.iat[0,5]/1000), (multiple_new.iat[24,6])/(multiple_new.iat[0,6]/1000)),
}
dic_lifetime_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[5,5], single_old.iat[5,6]),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[5,5], multiple_old.iat[5,6]),
    "single_new_price": linear_interpolate_2045(single_new.iat[5,5], single_new.iat[5,6]),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[5,5], multiple_new.iat[5,6]),
}

dic_marginal_cost_2045 = {
    "single_old_price": linear_interpolate_2045(single_old.iat[23,2] / 1000, single_old.iat[23,2] / 1000),
    "multiple_old_price": linear_interpolate_2045(multiple_old.iat[23,2] / 1000, multiple_old.iat[23,2] / 1000),
    "single_new_price": linear_interpolate_2045(single_new.iat[23,2] / 1000,single_new.iat[23,2] ),
    "multiple_new_price": linear_interpolate_2045(multiple_new.iat[23,2] / 1000, multiple_new.iat[23,2] / 1000),
}

dic_2020 = [dic_capacity_cost_overnight_2020, dic_effiency_2020, dic_fixom_cost_2020, dic_lifetime_2020, dic_marginal_cost_2020]
dic_2045 = [dic_capacity_cost_overnight_2045,dic_effiency_2045, dic_fixom_cost_2045, dic_lifetime_2045, dic_marginal_cost_2045]
val_2020 = []
val_2045 = []

# Berechnungen
for dic in dic_2020:
    val_2020.append(get_agg_price_2020(dic))

for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2020, val_2045)

```
