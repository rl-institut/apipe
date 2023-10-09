# 'Raw' Datasets 

------------------------------
## Technologiedaten

### Jahresvolllaststunden

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

### Leistungsdichte

Installierbare Leistung pro Fläche / spezifischer Flächenbedarf:

- Windenergie: 21 MW/km²
- PV-Freiflächenanlagen: 100 MW/km²
- PV-Aufdachanlagen: 140 MW/km²
- Solarthermie: ? MW/km²

Quelle: [PV- und Windflächenrechner](https://zenodo.org/record/6794558)

Datei: `technology_data.json` --> `power_density`

### Nennleistung Windenergieanlage

Als Zukunftsanlage für 2045 wird eine Enercon E126 6500 (6,5 MW) angenommen.
Diese wird für die Berechnung der Anlagenanzahl in den Ergebnissen
verwendet.

Datei: `technology_data.json` --> `nominal_power_per_unit`

### Batterien

- Kleinbatterien/Heimspeicher: Nennkapazität je installierter PV-Peakleistung
  und Speichernennleistung je installierter Speichernennkapazität aus
  [bnetza_mastr](../../digipipe/store/digipipe/store/raw/bnetza_mastr/dataset.md) und
  [HTW](https://solar.htw-berlin.de/wp-content/uploads/HTW-Stromspeicher-Inspektion-2023.pdf).
- Großbatterien: Speichernennleistung je installierter Speichernennkapazität
  aus [bnetza_mastr](../../digipipe/store/digipipe/store/raw/bnetza_mastr/dataset.md).

Datei: `technology_data.json` --> `batteries`

### Warmwasserspeicher

- Kleinwärmespeicher (dezentral): Speichernennleistung je installierter
  Speichernennkapazität aus
  [DEA](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage)
- Großwärmespeicher (Fernwärme): Speichernennleistung je installierter
  Speichernennkapazität aus
  [DEA](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage)

Datei: `technology_data.json` --> `hot_water_storages`

### Kosten und Wirkungsgrade

Datei: `raw_costs_efficiencies.csv`

##### Allgemein

Preise werden aus den Technologie Datenblättern der Danish Energy
Agency ([1], [2], [3], [4]) entnommen.
Abweichungen werden gesondert genannt.

alle Preise werden auf Euro im Jahr 2020 (dis-)kontiert und damit
inflationsbereinigt.

Für Quellen
[1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and),
[2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants),
[3](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage),
[4](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-renewable-fuels)
ist das meist die Umrechnung von 2015 zu 2020. Dafür folgende Formel verwendet:

```
P_(2020) = P_(2015)*f_(infl)
f_(infl) = (1+i_(2015))*(1+i_(2016))...*(1+i_(2019))
f_(infl) = 1,005 * 1,005 * 1.015 * 1,018 * 1,014 = 1,0582
```

[8](https://de.statista.com/themen/112/inflation/#topicOverview)

Werte für 2045 werden durch lineare Extrapolation ermittelt.

##### biogas_upgrading plant

Quelle: [4] "82 Biogas, upgrading"

Aufbereitung von Biogas zu Bio-SNG

##### biogas bpchp_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
"06 Gas engines, biogas"

Backpressure Combined heat and power (bpchp) modelliert BHKWs

thermal effiency = electrical_effiency / (c_b+c_v)  (
laut [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
S. 390)

##### biogas bpchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
"06 Gas engines, biogas"

Identische Werte zu biogas bpchp_central. Split fürs Energiesystem, aber
eingesetzte Technologie identisch.

##### biogas_plant

Quelle [4](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-renewable-fuels):
"81 Biogas Plant, Basic conf."

Stellt Biogas bereit, welches in KWK (biogas bpchp_central, biogas
bpchp_decentral) genutzt werden kann

##### boiler_central

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and):
"44 Natural Gas DH Only"

##### boiler_decentral

Quelle [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants):
"202 Gas boiler, ex single", "202 Gas boiler, ex apart", "202 Gas boiler, new
single",
"202 Gas boiler, new apart"

Es werden für jedes Szenario jeder Wert aus 4 Komponenten zusammengesetzt.

Diese sind die Kombinationen aus:

- Altbau-Neubau
- Einfamilienhaus-Mehrfamilienhaus

Diese Kompnonten werden durch Faktoren gewichtet zusammengefasst.

Für 2020:

- Verhältnis von Altbau-Neubau
  aus [7](https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/)
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Neubau
  aus [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb),
  verbaute Gasheizungen aggregiert
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Altbau wird als 0.7 / 0.3
  angenommen

Für 2045:

- Verhältnis von Altbau-Neubau
  aus [7](https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/)
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Neubau
  aus [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb),
  verbaute Gasheizungen in 2020
- Verhätnis von Einfamilienhaus-Mehrfamilienhaus im Altbau wird als 0.7 / 0.3
  angenommen

volle Berechnungen siehe "boiler_small_script.py" im Code Anhang

##### ch4 bpchp_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
"06 Gas engines, natural gas"

Backpressure Combined heat and power (bpchp) modelliert BHKWs

thermal effiency = electrical_effiency / (c_b+c_v)  (
laut [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
S. 390)

##### ch4 bpchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"06 Gas engines, natural gas"

Identische Werte zu ch4 bpchp_central. Split fürs Energiesystem, aber
eingesetzte Technologie identisch.

##### ch4 extchp_central

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"05 Gas turb. CC, steam extract., Large", [14] S. 20-21

##### ch4 extchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"05 Gas turb. CC, steam extract., Large"

[14] S. 20-21

Identisch wie ch4 extchp_central

##### gt

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"
04 Gas turb. simple cycle, L"

gas turbine, offener Prozess

##### heatpump_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"
40 Comp. hp, airsource 10 MW"

Wärmepumpentechnologie (Luft-Wasser-WP) aus Langfristigkeitsszenarien

##### heatpump_decentral

Quellen: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants
"
207 HP air-water,ex single", "207 HP air-water,ex apart", "207 HP air-water,new
single", "207 HP air-water,new apart", "
207 HP ground-water,ex single",  "207 HP ground-water,ex apart", "207 HP
ground-water,new single", "207 HP
ground-water,new apart",
[5], [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb)

Es werden für jedes Szenario jeder Wert aus 8 Komponenten zusammengesetzt.
Diese sind die Kombinationen aus:

- Sole-Umwelt
- Einfamilienhaus-Mehrfamilienhaus (fast alle WP in Einfamilienhäsuern!)
- Altbau-Neubau

Es wird das gemittelte Verhätnis Deutschlandweit der letzten 20 Jahre
angenommen (BBSR; Bundesamt für Bauwesen und
Raumordnung)

Für 2020 wurden Annahmen für das allgemeine Verhältnis zwischen den
Möglichkeiten angenommen:

- Sole-Umwelt sind die aggregierten Absatzzahlen aus [5]
- Einfamilienhaus-Mehrfamilienhaus
  aus [6](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb)
- Altbau-Neubau
  aus [7](https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/)

Mit diesen wird für 2045 wurden Annahmen für das allgemeine Verhältnis zwischen
den Möglichkeiten angenommen:

- Sole-Umwelt = 0.87/0.13 (Das sind die Absatzzahlen aus 2022 aus der
  Branchenstudie)
- Einfamilienhaus-Mehrfamilienhaus = 0.7 / 0.3 (Das ist eine freie Annahme, die
  eine fortschreitende Verbreitung in
  Mehrfamilienhäusern annimmt)
- Altbau-Neubau = 0.699 / 0.301 (das gemittelte Verhätnis Deutschlandweit der
  letzten 20 Jahre)

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

##### large_scale_battery

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"
180 Lithium Ion Battery", Cebulla [9] S. 181

storage_fixom_cost Berechnung aus UMAS/Oemof_B3 übernommen, ohne Quelle dieser
Berechnung gefunden zu haben.

storage_fixom_cost = 0,005 * storage_capacity_cost_overnight

Große Differenzen zwischen Windnode und UMAS, UMAS Methodik übernommen

##### pth_central

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"
41 Electric Boilers, small", "41 Electric Boilers, large"

Es wurde ein Mittelwert aus den Electric Biolers small und large gebildet, um
relevante Größen in ABW abzubilden.

##### pth_decentral

Quellen: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants): "
216 Electric heating,new single", "216 Electric heating,new apart"

Annahmen zu Gebäudebestand siehe heatpump_decentral, nur ohne Kombination mit
Altbau, da power to heat in Altbauten
vernachlässigbar selten (und wenn in anderen Technologien wie
Nachtspeicherheizungen) vorkommt.

Berechnungen siehe "pth_decentral_script" im Code Anhang

##### small_scale_battery

Quelle: [15](https://www.zhb-flensburg.de/fileadmin/content/spezial-einrichtungen/zhb/dokumente/dissertationen/fluri/fluri-2019-wirtschaftlichkeit-dez-stromspeicher.pdf), [17]
S. 3

- capacity_cost_overnight: [15] S. 41
- effiency, lost_rate, lifetime: [15] S.v91

##### storage heat_central

Quelle [3](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage): "
141
Large hot water tank"

- capacity_cost_overnight und fixom_cost ignoriert, da
  storage_capacity_cost_overnight, storage_fixom_cost einen Wert
  hat
- storage heat_decentral

Quelle [3](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage): "
141
Large hot water tank"

capacity_cost_overnight und fixom_cost ignoriert, da
storage_capacity_cost_overnight, storage_fixom_cost einen Wert hat

Große Differenzen zwischen UMAS und Windnode, UMAS Methodik übernommen

##### hydro ror

Quellen: [16]

- fixom_cost: S. 78
- capacity_cost_overnight: S.75
- lifetime: S. 72

##### lignite oven

Quellen: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)  "
206 Wood stove, single, ex tank"

Der Kohleofen ist eine Komponente, die für die Abbildung des Ist-Zusandes
relevant ist.
Die Kohleheizung wird durch gesetzliche Regulierung nicht mehr neu verbaut
werden können, wodurch die Komponente für die
Optimierung nicht relevant ist.
Auch die Datenlage für die Kohleheizung sehr schlecht ist, die Daten werden
daher approximiert.

Keine direkten Werte vorhanden, daher Modellierung anhand der wood stove Werte

efficiency:

Differenz der Energie zwischen Holz und Kohle liegt im Heizwert des Brennstoffs.
Daher wird die Effizienz der wood stove
mit Faktor des Verhältnisses der Heizwerte multipliziert.
Daten für Heizwerte von
BMWK [11](https://www.bmwk.de/Redaktion/DE/Artikel/Energie/energiedaten-gesamtausgabe.html)
und [12](https://books.google.de/books?id=n0fVYjrHAlwC&pg=PA58#v=onepage&q&f=false)
ergibt einen Faktor von 4/3

fixom_cost:

Bestehen großteils aus Brennstoffkosten. Änderung zu wood stove besteht aus
Heizwert (gewonnene Energie pro kg) und
Preisdiff pro Kilogramm

Preise aus brikett-rekord.com [13]

lifetime:

identisch wie wood stove

marginal-cost: identisch wie wood stove

Aus den Annahmen folgt, dass die Investkosten ignoriert werden können.

##### pv_ground

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "
22 Utility-scale PV", Vergleich [10]

marginal_cost = 0, da in Quellen nicht vorhanden

Kosten
aus [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
im Bereich von [10]

##### pv_rooftop

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "
22 PV commercial&industrial rooftop", "22 PV residential", Vergleich [10]

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

Kosten
aus [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
im Bereich von [10]

##### thermalcollector_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and
"
46 Solar District Heating"

##### thermalcollector_decentral

Quelle: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants
"
215 Solar heating,ex single", "215 Solar heating,ex apart", "215 Solar
heating,new single", "215 Solar heating,new
apart"

Annahmen zu Gebäudebestand siehe heatpump_decentral.

Berechnungen siehe "thermalcollector_decentral_script" im Code Anhang

##### wind onshore

Quelle [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and): "
20 Onshore turbines",
Vergleich [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)

EE Kosten durchweg kleiner als in Windnode in 2020

Windnode bezieht sich auf Frauenhofer ISE aus 2018, Vorgängerstudie
zu [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)

Frauenhofer (S.

11) [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)
    CAPEX-Range höher als
    DEA [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
    in 2020

1400000-2000000 [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)
zu
1190000 [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)
€/MW

keine Aussagen in
Frauenhofer [10](https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html)
über 2045

wir wählen DEA als Quelle für die Vergleichbarkeit, da Vergleichbarkeit in der
Optimierung der Modellierung Vorrang hat

##### wood extchp_central

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)  "
09a Wood Chips, Medium"

[14] S. 20-21

##### wood extchp_decentral

Quelle: [1](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and)  "
09a Wood Chips, Medium"

[14] S. 20-21

identisch zu wood extchp_central

##### wood oven

Quelle: [2](https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants), "
204 Biomass auto,ex single", "204 Biomass auto,new single", "204 Biomass auto,ex
apart", "204 Biomass auto,new apart"

Annahmen zu Gebäudebestand siehe heatpump_decentral.

Berechnungen siehe "wood_oven_script" im Code Anhang

##### Quellen

[1] Danish Energy Agency (2016): "Technology Data - Energy Plants for
Electricity and District heating generation",
Version 13,
von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-generation-electricity-and

[2] Danish Energy Agency (2016): "Technology Data for heating installations",
Version 4,
von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-individual-heating-plants

[3] Danish Energy Agency (2018): "Technology Data – Energy storage", Version 7,
von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-energy-storage

[4] Danish Energy Agency (2017): "Technology Data – Renewable fuels", Versoin 9,
von https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-renewable-fuels

[5] Karl-Heinz Backhaus (Vaillant), Dr. Hendrik Ehrhardt (Stiebel Eltron), Sven
Kersten (NIBE), Steffen
Moser (EnBW), Frank Richert (Wolf), Ingo Rieger (Bosch), Egbert Tippelt (
Viessmann), André Jacob
(BWP), Johanna Otting (BWP), Björn Schreinermacher (BWP)(2023): "Branchenstudie
2023: Marktentwicklung – Prognose
–Handlungsempfehlungen", Bundesverband Wärmepumpe (BWP) e. V.

[6] Statistisches Landesamt Sachsen-Anhalt: "GENESIS-Tabelle: 31121-0006,
Statistik der Baufertigstellungen",
von https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=31121-0006&bypass=true&levelindex=0&levelid=1682324189765#abreadcrumb,
Stand: 11.04.2023

[7] Statista Research Department(2021): "Struktur des Wohnungsbaus nach Neubau
und Sanierung in Deutschland in den
Jahren 2001 bis 2020",
von https://de.statista.com/statistik/daten/studie/202207/umfrage/struktur-des-wohnungsbaus-nach-art-der-bauleistung-in-deutschland/,
Stand: 03.04.2023 12:26:20

[8] Statista: "Daten und Fakten zur Inflation und den Verbraucherpreisen" ,
von https://de.statista.com/themen/112/inflation/#topicOverview , Stand:
29.03.2023

[9] Cebulla, Felix (2017): "Storage demand in highly renewable energy scenarios
for Europe", OPUS - Online Publikationen
der Universität Stuttgart, von https://elib.uni-stuttgart.de/handle/11682/9778

[10] Frauenhofer ISE (2019): "Stromgestehungskosten erneuerbare Energien",
von https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html

[11] BMWK (2021): "Energiedaten"
von https://www.bmwk.de/Redaktion/DE/Artikel/Energie/energiedaten-gesamtausgabe.html

[12] Michael Herrmann, Jürgen Weber: Öfen und Kamine: Raumheizungen fachgerecht
planen und bauen. Beuth Verlag, 201,
von  https://books.google.de/books?id=n0fVYjrHAlwC&pg=PA58#v=onepage&q&f=false

[13] www.brikett-rekord.com: "Energiekostenvergleich",
von https://www.brikett-rekord.com/de/heizwertvergleich-rekord-briketts.html,
letzter Abruf 8.5.2023

[14] WindNode: Modell, Methodik, Daten, ABW;
von: RLI
letzer Abruf 8.8.2023

[15] Fluri, Verena: "Wirtschaftlichkeit von zukunftsfähigen Geschäftsmodellen
dezentraler Stromspeicher"
von https://www.zhb-flensburg.de/fileadmin/content/spezial-einrichtungen/zhb/dokumente/dissertationen/fluri/fluri-2019-wirtschaftlichkeit-dez-stromspeicher.pdf,
letzter Abruf 8.8.2023

[16] Schröder, Andreas; Kunz, Friedrich; Meiss, Jan; Mendelevitch, Roman;
Hirschhausen, Christian von: "Current and Prospective Costs of Electricity
Generation until 2050"
von https://www.diw.de/documents/publikationen/73/diw_01.c.424566.de/diw_datadoc_2013-068.pdf,
letzter Abruf 8.8.2023

[17] Prüggler, Wolfgang (2019): "HEIMSPEICHERSYSTEME UND ZENTRALE
BATTERIESPEICHER – KRITISCHE FAKTOREN DER WIRTSCHAFTLICHKEIT"
von https://ens.dk/sites/ens.dk/files/Analyser/technology_data_catalogue_for_energy_storage.pdf,
letzter Abruf 8.8.2023

#### Code Anhang

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


## Daten aus DEA:
## einlesen von Daten
data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

##datensheets
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

## Berechnungen
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


## Daten aus DEA:
## einlesen von Daten
data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

##datensheets
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

## Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2045)
```

pv_rooftop_script.py:

```
import pandas as pd
import geopandas as gpd
import os.path

##trennt residential and industrial rooftop PV nach Nennleistung
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

## Daten aus DEA:
## einlesen von Daten
data_folder_sheets = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data_sheets = os.path.join(data_folder_sheets, "technology_data_for_el_and_dh.xlsx")

##datensheets
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

## Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic, proportion))

print(dic_capacity_cost_overnight_2045, dic_fixom_cost_2045, dic_lifetime_2045)
print(proportion[0][0])
print(val_2045)
```

Pth_decentral_sc0irpt.py

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


## Daten aus DEA:
## einlesen von Daten
data_folder = os.path.join("/YOUR/DATA/ROAD/TAKE/ME/HOME/TO/THE/PLACE")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

##datensheets
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
## Berechnungen
for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2045)
```

heatpump_small_script.py:

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


## Daten aus DEA:
## einlesen von Daten
data_folder = os.path.join("/home/local/RL-INSTITUT/aaron.schilling/Dokumente/Projekte/Digipipe")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

##datensheets
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

## Berechnungen
for dic in dic_2020:
    val_2020.append(get_agg_price_2020(dic)[0])

for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic)[0])

print(val_2020, val_2045)
print(faktoren_2020, faktoren_2045)
```

boiler_small_script.py:

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


## Daten aus DEA:
## einlesen von Daten
data_folder = os.path.join("/ROAD/TO/DATA")
data = os.path.join(data_folder, "technology_data_heating_installations_-_0003.xlsx")

##datensheets
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

## Berechnungen
for dic in dic_2020:
    val_2020.append(get_agg_price_2020(dic))

for dic in dic_2045:
    val_2045.append(get_agg_price_2045(dic))

print(val_2020, val_2045)

```

**Dataset: `raw/technology_data`**

??? metadata "Metadata"
    ```json
    {
        "name": "technology_data",
        "title": "Technologiedaten",
        "id": "technology_data",
        "description": "Jahresvollaststunden, Leistungsdichte, Nennleistung, Kosten und Effizienzen von energieumwandlungs Technologien",
        "language": [
            "de-DE",
            "en-GB"
        ],
        "subject": [],
        "keywords": [
            "technologiedaten",
            "Jahresvollaststunden",
            "Nennleistung",
            "Kosten",
            "Effizienzen"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Europe",
            "extent": "Europe",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": null,
            "timeseries": null
        },
        "sources": [
            {
                "title": "F\u00f6deral Erneuerbar",
                "description": "Jahresvollaststunden, Leistungsdichte, Nennleistung, Kosten und Effizienzen von energieumwandlungs Technologien",
                "path": "https://www.foederal-erneuerbar.de",
                "licenses": null
            },
            {
                "title": "PV- und Windfl\u00e4chenrechner",
                "description": "Der Photovoltaik- und Windfl\u00e4chenrechner - Methoden und Daten",
                "path": "https://zenodo.org/record/6794558",
                "licenses": null
            },
            {
                "title": "Ariadne Szenarienreport",
                "description": "Ariadne Szenarienreport",
                "path": "https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf",
                "licenses": null
            },
            {
                "title": "Technologiedaten",
                "description": "Kosten und Effizienzen von Energieumwandlungtechnologien",
                "path": "https://ens.dk/en/our-services/projections-and-models/technology-data",
                "licenses": null
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": [],
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Lokale Verwaltungseinheiten

Lokale Verwaltungseinheiten (LAUs) von Eurostat, mit NUTS kompatibel. Diese LAU
sind die Bausteine der NUTS und umfassen die Gemeinden und Kommunen der
Europäischen Union.

**Dataset: `raw/eurostat_lau`**

??? metadata "Metadata"
    ```json
    {
        "name": "eurostat_lau",
        "title": "Lokale Verwaltungseinheiten",
        "id": "eurostat_lau",
        "description": "Lokale Verwaltungseinheiten (LAUs) von Eurostat, mit NUTS kompatibel. Diese LAU sind die Bausteine der NUTS und umfassen die Gemeinden und Kommunen der Europ\u00e4ischen Union",
        "language": [
            "en-GB"
        ],
        "subject": null,
        "keywords": [
            "Verwaltungseinheiten",
            "NUTS",
            "LAU"
        ],
        "publicationDate": "2022-12-15",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": "NUTS-3"
        },
        "temporal": {
            "referenceDate": "2022-06-30",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Lokale Verwaltungseinheiten",
                "description": "Lokale Verwaltungseinheiten (LAUs) von Eurostat, mit NUTS kompatibel. Diese LAU sind die Bausteine der NUTS und umfassen die Gemeinden und Kommunen der Europ\u00e4ischen Union",
                "path": "https://ec.europa.eu/eurostat/de/web/nuts/local-administrative-units",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 eurostat, 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 eurostat, 2023"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-25",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Bevölkerungsprognose Sachsen-Anhalt

Bevölkerungsprognose je Gemeinde bis 2035 des Statistischen Landesamtes
Sachsen-Anhalt. Stand: 2021

**Dataset: `raw/stala_st_pop_prog`**

??? metadata "Metadata"
    ```json
    {
        "name": "stala_st_pop_prog",
        "title": "Regionalisierte Bev\u00f6lkerungsprognose",
        "id": "stala_st_pop_prog",
        "description": "Prognostizierter Bev\u00f6lkerungsstand in den Gemeinden, kreisfreien St\u00e4dten und Landkreisen nach Prognosejahr und Geschlecht",
        "language": [
            "de-DE"
        ],
        "subject": [],
        "keywords": [
            "Bev\u00f6lkerungsprognose",
            "population"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Sachsen-Anhalt",
            "extent": "Sachsen-Anhalt",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": null,
            "timeseries": [
                {
                    "start": "2019",
                    "end": "2035",
                    "resolution": "1 year",
                    "alignment": null,
                    "aggregationType": "sum"
                }
            ]
        },
        "sources": [
            {
                "title": "1_Internettabelle_7RBP_nach_Prognosejahr_Geschlecht_alle_Ebenen",
                "description": "Prognostizierter Bev\u00f6lkerungsstand in den Gemeinden, kreisfreien St\u00e4dten und Landkreisen nach Prognosejahr und Geschlecht",
                "path": "statistik.sachsen-anhalt.de/themen/bevoelkerung-mikrozensus-freiwillige-haushaltserhebungen/bevoelkerung/bevoelkerungsprognose-und-haushalteprognose/#c312231",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 Version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 2023 Landesportal Sachsen-Anhalt "
                    }
                ]
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": [],
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft.

Dafür wurden auf Basis des
[Digitalen Oberflächenmodells (DOM2)](https://www.lvermgeo.sachsen-anhalt.de/de/dom2-landesweit.html)
Schattenberechnungen durchgeführt. Anhand des
[LoD2 3D-Gebäudemodells](https://www.lvermgeo.sachsen-anhalt.de/de/download_lod2.html)
wurden für verschiedene Dachausrichtungen (nord, ost, süd, west, flach) die
installierbare Leistung bestimmt und mittels der Globalstrahlung und typischer
technischer Parameter für jedes Gebäude und jede Dachflächenorientierung
potenzielle Erträge berechnet.

Quellen:

- [Hauptseite](https://www.planungsregion-abw.de/geodaten/)
- [Geodaten](https://gis-entwicklung2.planungsregion-abw.de/geoserver/wfs?SERVICE=WFS&REQUEST=GetCapabilities)
- [Anwendung](https://ris.planungsregion-abw.de/mapbender/application/pv_dachflaechenpot_rpg_abw)

**Dataset: `raw/rpg_abw_pv_roof_potential`**

??? metadata "Metadata"
    ```json
    {
        "name": "rpg_abw_pv_roof_potential",
        "title": "Dachfl\u00e4chenpotenzial PV-Aufdachanlagen in ABW",
        "id": "rpg_abw_pv_roof_potential",
        "description": "Absch\u00e4tzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft.",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "PV",
            "Photovoltaic",
            "Fl\u00e4chenpotential",
            "Aufdachanlagen"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "ABW",
            "extent": "ABW",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-06-10",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Dachfl\u00e4chenpotenzial PV-Aufdachanlagen in ABW",
                "description": "Absch\u00e4tzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft.",
                "path": "https://www.planungsregion-abw.de/geodaten/",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Regionale Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg / 2022"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Regionale Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg / 2022"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Installierte Leistungen von Biomasse-Konversionstechnologien

Die installierten Leistungen in MW wird im Szenario 80 % Transformationspfad
und 2,6 Mio. ha Anbauflächen im Jahr 2020 und 2050 der Tabelle 13 im
Dokument
["Technoökonomische Analyse und Transformationspfade des energetischen Biomassepotentials (TATBIO)"](../../digipipe/store/raw/dbfz_biomass_heat_capacities/metadata.json)
für die folgenden Konversionsanlagen von Biomasse entnommen:

- Biomethan-Blockheizkraftwerk
- Holzhackschnitzelkessel Sektor Industrie
- Pelletkessel Sektor GHD
- Holzhackschnitzelkessel Sektor GHD
- Scheitholzvergaserkessel
- Pelletkessel Sektor Gebäude
- Biogasanlage + Blockheizkraftwerk
- Biomethan Gas- und Dampfkombikraftwerk
- Klärschlammfaulung + Blockheizkraftwerk
- Papier-Zellstoff-KWK
- Holzvergaser + Blockheizkraftwerk
- Mikro-Holzgas-Blockheizkraftwerk

Die Konversionstechnologien sind in der Spalte "technology" gelistet, während
sich ihre installierten Leistungen für die beiden Projektionsjahre in den
Spalten "capacity_[MW]_2020" und "capacity_[MW]_2050" befinden.

In den Spalten "decentral" und "central" wird mit "x" angegeben, ob jeweils ein
dezentraler und zentraler Einsatz der Konversionsanlage Stand der Technik ist.

In der Spalte "carrier" wird analog zur Konvention der Namensgebung im
Energiesystem (siehe [esys.md](../../digipipe/store/../../docs/sections/esys.md)) der
jeweilige in die Konversionsanlage eintretende Energieträger notiert.
Diese werden Abbildung 3 des Dokuments entommen. Der Energieträger Schwarzlauge
wird vereinfachend dem Energieträger feste Biomasse bzw. Holz zugeordnet.
Klärgas und Holzgas werden vereinfachend Biogas zugeordnet.

In der Spalte "tech" findet die Zuordnung zu der Technologie anhand der im
Energiesystem verwendeten Komponenten (siehe
[esys.md](../../digipipe/store/../../docs/sections/esys.md)) statt.

**Dataset: `raw/dbfz_biomass_heat_capacities`**

??? metadata "Metadata"
    ```json
    {
        "name": "dbfz_biomass_heat_capacities",
        "title": "Techno\u00f6konomische Analyse und Transformationspfade des energetischen Biomassepotentials (TATBIO)",
        "id": "dbfz_biomass_heat_capacities",
        "description": "Installierte Leistungen von Biomasse-Konversionstechnologien. Die installierten Leistungen in MW wird im Szenario 80 % Transformationspfad",
        "language": [
            "de-DE"
        ],
        "subject": [],
        "keywords": [
            "Biomasse",
            "Biomassepotential",
            "Analyse",
            "Transformationspfade",
            "TATBIO"
        ],
        "publicationDate": "2019-05-08",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": "2019-04-30",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Techno\u00f6konomische Analyse und Transformationspfade des energetischen Biomassepotentials (TATBIO)",
                "description": "Installierte Leistungen von Biomasse-Konversionstechnologien. Die installierten Leistungen in MW wird im Szenario 80 % Transformationspfad",
                "path": "https://www.ufz.de/export/data/2/231891_technooekonomische-analyse-und-transformationspfade-des-energetischen-biomassepotentials(1).pdf",
                "licenses": null
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": [],
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Temperatur

Stündliche Mittelwerte der Luft- und Erdbodentemperatur des Deutschen
Wetterdienstes
([Climate Data Center](https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/))
für das Jahr 2011 je Gemeinde in der Region ABW, vorverarbeitet im Projekt
[WindNODE](https://windnode-abw.readthedocs.io/en/latest/energy_system_model.html#energy-demand-today).

Werte:

- `temp_amb`: Lufttemperatur in 2 m Höhe
- `temp_soil`: Erdbodentemperatur in 1 m Tiefe

Verwendete Stationen:

- Wittenberg
- Köthen
- Jessnitz
- Seehausen
- Holzdorf

Die Zuordnung der Stationsmesswerte zu Gemeinden erfolgte über die jeweils
nächstgelegene Wetterstation.

**Dataset: `raw/dwd_temperature`**

??? metadata "Metadata"
    ```json
    {
        "name": "dwd_temperature",
        "title": "temperatur_2011",
        "id": "dwd_temperature",
        "description": "St\u00fcndliche Mittelwerte der Luft- und Erdbodentemperatur des Deutschen Wetterdienstes (Climate Data Center) f\u00fcr das Jahr 2011 je Gemeinde in der Region ABW, vorverarbeitet im Projekt WindNODE.",
        "language": [
            "en-GB"
        ],
        "subject": null,
        "keywords": [
            "Wetter",
            "Temperatur"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Anhalt-Bitterfeld-Wittenberg",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2023-05-23",
            "timeseries": [
                {
                    "start": "2011-01-01T00:00+01",
                    "end": "2011-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": "left",
                    "aggregationType": null
                }
            ]
        },
        "sources": [
            {
                "title": "temperatur_2011",
                "description": "St\u00fcndliche Mittelwerte der Luft- und Erdbodentemperatur des Deutschen Wetterdienstes (Climate Data Center) f\u00fcr das Jahr 2011 je Gemeinde in der Region ABW, vorverarbeitet im Projekt WindNODE",
                "path": "https://www.dwd.de/DE/leistungen/cdc/climate-data-center.html",
                "licenses": [
                    {
                        "name": "GeoNutzV",
                        "title": "Verordnung zur Festlegung der Nutzungsbestimmungen f\u00fcr die Bereitstellung von Geodaten des Bundes",
                        "path": "https://www.gesetze-im-internet.de/geonutzv/GeoNutzV.pdf",
                        "instruction": "Alle frei zug\u00e4nglichen Geodaten und Geodatendienste d\u00fcrfen entsprechend der Verordnung zur Festlegung der Nutzungsbestimmungen f\u00fcr die Bereitstellung von Geodaten des Bundes (GeoNutzV) unter Beigabe eines Quellenvermerks ohne Einschr\u00e4nkungen weiterverwendet werden.",
                        "attribution": " \u00a9 Deutscher Wetterdienst"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "GeoNutzV",
                "title": "Verordnung zur Festlegung der Nutzungsbestimmungen f\u00fcr die Bereitstellung von Geodaten des Bundes",
                "path": "https://www.gesetze-im-internet.de/geonutzv/GeoNutzV.pdf",
                "instruction": "Alle frei zug\u00e4nglichen Geodaten und Geodatendienste d\u00fcrfen entsprechend der Verordnung zur Festlegung der Nutzungsbestimmungen f\u00fcr die Bereitstellung von Geodaten des Bundes (GeoNutzV) unter Beigabe eines Quellenvermerks ohne Einschr\u00e4nkungen weiterverwendet werden.",
                "attribution": " \u00a9 Deutscher Wetterdienst"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "aaron.schilling@rl-institut.de",
                "date": "2023-08-25",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Emissionen

Emissionen für die Jahre 1990 und 2019 für Sachsen-Anhalt und disaggregiert für
die Region Anhalt-Bitterfeld-Wittenberg (ABW). Die Grundlage hierfür ist der
[THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf)
Sachsen-Anhalt (ST).

Datei: `emissions.csv`, Felder:

- `sector`: Sektor
- `cat`: Kategorie ("*" = alle)
- `subcat`: Unterkategorie ("*" = alle)
- `name`: Bezeichner
- `st`: Emissionen Sachsen-Anhalt in kt CO2-Äquivalent
- `abw`: Emissionen Region ABW in kt CO2-Äquivalent

`sector`, `cat` und `subcat` folgen der Nomenklatur des Common Reporting Formats
(CRF) nach [KSG Anlage 1](https://www.gesetze-im-internet.de/ksg/anlage_1.html).
[Grafik hierzu](https://expertenrat-klima.de/content/uploads/2023/05/ERK2023_Pruefbericht-Emissionsdaten-des-Jahres-2022.pdf)
(Abb. 2 auf S. 30).

### Disaggregation

Anhand unterschiedlicher Kriterien und Datenquellen wurde näherungsweise von den
vorliegenden Emissionen für Sachsen-Anhalt für 1990 und 2019 auf die Region ABW
disaggregiert. Je Sektor sind hier die gewählten
**energiebestimmenden Größen (EnbG)** angegeben, sowie die Herangehensweise zur
jeweiligen Berechnung.

#### Sektor Energiewirtschaft (CRF 1.A.1 + 1.B)

Aus der Liste der
[Emissionshandelspflichtigen Anlagen](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3)
wurden jene Daten zu Anlagen extrahiert, welche sich in Sachsen-Anhalt befinden
und als Bezeichnung "Energieumwandlung >= 50 MW FWL" oder "Energieumwandlung
20–50 MW FWL" (Haupttätigkeit nach TEHG) aufweisen.
Die Summe der angegebenen Emissionen (t CO2 Äq) jener Anlagen, welche in der
Region ABW liegen, wurde in Relation zu der Summe der Emissionen aus den Anlagen
in Gesamt-ST gesetzt. Dieser Anteil wurde auf die im THG-Bericht angegebene
Emissionsmenge im Sektor "Energiewirtschaft (1.A.1)" sowie "Prozessemissionen
(1.B)" angelegt und so für ABW näherungsweise disaggregiert.

Hinweise:

- Aufgrund mangelnder Daten wurde für das Jahr 1990 auf die neuesten verfügbaren
  Daten (2005-2007) aus der Anlagenliste zurückgegriffen.
- Energiewirtschaftlich relevante Anlagen unter 20 MW FWL sind in der
  Anlagenliste nicht erfasst und konnten somit nicht berücksichtigt werden.

Quellen:

- [Emissionshandelspflichtige Anlagen in Deutschland 2020 (Stand 03.05.2021)](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3)
- [Treibhausgasemissionen in Sachsen-Anhalt 2018 (Stand 12.05.2021)](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/THG_Bericht_2018.pdf)

##### CRF 1.A.1

Energiewirtschaft (Umwandlungsbereich): umfasst die öffentliche Elektrizitäts-
und Wärmeversorgung sowie Raffinerien.

EnbG: Emissionen aus europäischem Emissionshandel

##### CRF 1.B

Diffuse Emissionen aus Brennstoffen: Diese Kategorie beinhaltet flüchtige
Emissionen aus der Gewinnung, Verarbeitung und Verteilung von Brennstoffen. Die
wichtigsten Quellen sind die Verteilung von Erdgas, aber auch Emissionen aus
Förderung und Abfackelung, die Extraktion und Umwandlung von Braunkohle,
Emissionen aus der Raffination von Erdöl sowie Emissionen aus der Lagerung und
Verteilung von Mineralölprodukten.

EnbG: Emissionen aus europäischem Emissionshandel

#### Sektor Industrie (CRF 1.A.2)

Dieser Sektor umfasst sämtliche energiebedingten Emissionen durch verarbeitendes
Gewerbe.

Zur Disaggregierung wurde der Energieverbrauch der Industriebetriebe in ABW mit
dem Gesamtenergieverbrauch aller Industriebetriebe in Sachsen-Anhalt in Relation
gesetzt. Dabei wurde eine Differenzierung hinsichtlich der
Energieträgerzusammensetzung von ABW im Vergleich zu ST durchgeführt und anhand
von Emissionsfaktoren berechnet.

EnbG: Energieverbrauch nach Energieträgern

Quellen:

- [Energieverbrauch der Industriebetriebe in Sachsen-Anhalt nach ausgewählten Energieträgern und Kreisen](https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Energie/Tabellen/Energieverwendung/Energieverbrauch_nach_Kreisen_ab_dem_Jahr_2010.xlsx)
- [Emissionsfaktor für Stromerzeugung (UBA)](https://www.umweltbundesamt.de/sites/default/files/medien/479/bilder/dateien/entwicklung_der_spezifischen_emissionen_des_deutschen_strommix_1990-2020_und_erste_schaetzungen_2021.pdf)
- [BISKO Bilanzierungs-Systematik Kommunal (Aktualisierung 11/2019)](https://www.ifeu.de/fileadmin/uploads/BISKO_Methodenpapier_kurz_ifeu_Nov19.pdf)

#### Sektor Prozessemissionen (CRF 2)

Dieser Sektor umfasst sämtliche Emissionen, welche durch Industrieprozesse
anfallen. Dies sind Emissionen aus: Herstellung mineralischer Produkte,
chemischer Industrie, Herstellung von Metallen, übrigen Prozessen und
Produktverwendungen (CRF 2.A-H).
Zur Disaggregierung wurde erneut die
[Liste der Emissionshandelspflichtigen Anlagen](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3)
herangezogen. Anders als im Sektor Energiewirtschaft (s.o.) wurde jedoch der
Anteil aller Anlagen, welche nicht der Energiewirtschaft zugerechnet werden, zur
Bestimmung des Anteils von ABW an ST gewählt.

EnbG: Emissionen aus europäischem Emissionshandel

#### Sektor Verkehr (CRF 1.A.3)

Dieser Sektor umfasst Emissionen aus dem Straßenverkehr, dem zivilen
Luftverkehr, aus dem Schiffsverkehr, verbrennungsbedingte Emissionen aus dem
Schienenverkehr sowie Emissionen des übrigen Verkehrs und weitere Quellen zur
Bereitstellung der im Verkehr verbrauchten Energie. Die Verbrennung von
Mineralölprodukten im Straßenverkehr spielt die größte Rolle und macht weit über
90 % der sektoralen Emissionen aus. Daher wird zur Disaggreagation der
motorisierte Straßenverkehr über zugelassene Kraftfahrzeuge mit
durchschnittlichen Fahrleistungen und spezifischer Emissionen pro Kilometer und
Fahrzeugklasse herangezogen.

Hierfür wird zunächst aus
[Verkehr in Kilometern (VK) ZeitreiheJahre 2014 - 2022](https://www.kba.de/DE/Statistik/Kraftverkehr/VerkehrKilometer/vk_inlaenderfahrleistung/vk_inlaenderfahrleistung_node.html;jsessionid=DD419FD0604C0BCC72A9E4533BB0319F.live21324)
und
[Umweltfreundlich mobil! Ein ökologischer Verkehrsartenvergleich für den Personen- und Güterverkehr in Deutschland)](https://www.umweltbundesamt.de/sites/default/files/medien/5750/publikationen/2021_fb_umweltfreundlich_mobil_bf.pdf)
ein durchschnittlicher Emissionswert pro Jahr und Fahrzeugklasse ermittelt.
Dieser wird mit den zugelassenen Fahrzeugen der entsprechenden Fahrzeugklassen
aus
[Kraftfahrzeugbestand nach Kraftfahrzeugarten - Stichtag 01.01. - regionale Tiefe: Kreise und krfr. Städte (bis 01.01.2019)](https://www-genesis.destatis.de/genesis//online?operation=table&code=46251-0001&bypass=true&levelindex=0&levelid=1691405772899#abreadcrumb)
einerseits für ganz Sachsen-Anhalt und andererseits ABW multipliziert. Daraus
wird ein Verhältnis der Verkehrsemissionen in ABW zu ST gewonnen.

Hinweise:

- Die Datenlage für die zugelassenen Fahrzeuge, gefahrenen Kilometer und
  Emissionen pro km sind nicht spezifisch für 1990 sondern nur für einzelne
  Jahre der frühen 1990er verfügbar. Daher ist der Emissionswert für 1990 mit
  einer höheren Unsicherheit behaftet.

EnbG:

- Zugelassene Kraftfahrzeuge
- Durchschnittliche Fahrleistung und spez. CO2 Emission pro km und
  Fahrzeugklasse

Quellen:

- [Kraftfahrzeugbestand nach Kraftfahrzeugarten - Stichtag 01.01. - regionale Tiefe: Kreise und krfr. Städte (bis 01.01.2019)](https://www-genesis.destatis.de/genesis//online?operation=table&code=46251-0001&bypass=true&levelindex=0&levelid=1691405772899#abreadcrumb)
- [Umweltfreundlich mobil! Ein ökologischer Verkehrsartenvergleich für den Personen- und Güterverkehr in Deutschland)](https://www.umweltbundesamt.de/sites/default/files/medien/5750/publikationen/2021_fb_umweltfreundlich_mobil_bf.pdf)
- [Verkehr in Kilometern (VK) ZeitreiheJahre 2014 - 2022](https://www.kba.de/DE/Statistik/Kraftverkehr/VerkehrKilometer/vk_inlaenderfahrleistung/vk_inlaenderfahrleistung_node.html;jsessionid=DD419FD0604C0BCC72A9E4533BB0319F.live21324)

#### Sektor Sonstige Energie (insbes. Gebäude) (CRF 1.A.4 + 1.A.5)

Dieser Sektor umfasst den durch Energieumwandlung nicht bereits abgedeckten
Energiebedarf. Das sind vor allem kleine Einzelfeuerungsanlagen bis hin zu
immissionsschutzrechtlich genehmigungsbedürftigen Anlagen mit einer
Nennwärmeleistung von mehreren Megawatt. Zur Disaggreagtion wurde daher der
Wärmebedarf von ABW im Verhältnis zum Wärmebedarf von gesamt Sachsen Anhalt
gewählt. Der Wärmevedarf umfasst Raumwärme, Warmwasser sowie Kochen und wird aus
Daten aus dem Pipeline-Datensatz
[demand_heat_region](../../digipipe/store/datasets/demand_heat_region/dataset.md) generiert.

Ergebnis: 17,46 % des Bedarfs in Sachsen-Anhalt entfällt auf ABW.

Code
```
## Sektor HH
heat_hh_dist_states = gpd.read_file("demand_heat_zonal_stats-res-bkg_vg250_federal_states.gpkg")
heat_hh_demand_st = float(heat_hh_dist_states.loc[heat_hh_dist_states.nuts == "DEE"].heat_demand)
heat_hh_demand_abw = gpd.read_file("demand_heat_zonal_stats-res-bkg_vg250_muns_region.gpkg").heat_demand.sum()

## Sektor GHD
heat_cts_dist_states = gpd.read_file("demand_heat_zonal_stats-ser-bkg_vg250_federal_states.gpkg")
heat_cts_demand_st = float(heat_cts_dist_states.loc[heat_cts_dist_states.nuts == "DEE"].heat_demand)
heat_cts_demand_abw = gpd.read_file("demand_heat_zonal_stats-ser-bkg_vg250_muns_region.gpkg").heat_demand.sum()

## Anteil ABW an ST
heat_share = (heat_hh_demand_abw + heat_cts_demand_abw) / (heat_hh_demand_st + heat_cts_demand_st)
```

EnbG: Wärmebedarf aus Energiesystem

#### Sektor Landwirtschaft (CRF 3)

Der Sektor umfasst Emissionen aus der Viehwirtschaft und der Bewirtschaftung von
Böden. Daher werden zunächst die Emissionsunterkategorien 3.A-J der
Viehwirtschaft oder der Bewirtschaftung von Böden zugeordnet. Anschließend
werden diese getrennt nach den Viehbeständen bzw. der landwirtschaftlich
genutzen Fläche disaggreiert.

##### CRF 3.A - Landwirtschaft – Fermentation

Emissionen durch Fermentation (CRF 3.A) entstehen vorrangig durch
Verdauungsprozesse in der Viehwirtschaft. Deswegen wird der Anteil ABWs an
diesen Emissionen durch die Viehbestände abgeschätzt.

Hinweis:

- Die Viehbestände für 1990 sind nicht bekannt, es wird stattdessen auf die
  Viehbestände von 1996 zurückggegriffen.

EnbG: Viehbestände

Quelle:

- [Viehbestand der landwirtschaftlichen Betriebe in Großvieheinheiten (GV) nach Jahren und Kreisen)](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/land-und-forstwirtschaft-fischerei/tabellen-viehwirtschaft-und-tierische-erzeugnisse#c234218)

###### CRF 3.B-J

In den Unterkategorien 3.C-J ist eine Proportionalität der Emissionen und der
landwirtschafltich genutzen Fläche zu erwarten. Unterkategorie 2.B
"Wirtschaftsdüngerausbringung (ohne Gärreste)" ist allerdings ein Grenzfall, da
er aus Abfällen der Tierhaltung produziert wird und bereits hierbei
Treibhausgase entstehen, diese aber nicht vor Ort eingesetzt werden müssen,
sondern auf beliebigen landwirtschafltichen Flächen eingesetzt werden kann.
Daher wird hier auch diese Unterkategorie der Landnutzung zugeordnet.

Hinweis:

- die Flächenntuzungsdaten gehen nicht bis 1990 zurück, ändern sich über die
  Jahre aber nur marginal, sodass hier nur von geringen Abweichungen auszugehen
  ist.

EnbG: Landwirtschaftlich genutzte Fläche

Quelle:

- [Flaeche_nach_Kultuarten_nach_Jahren_und_Kreisen](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/land-und-forstwirtschaft-fischerei/tabellen-bodennutzung-und-anbau)

#### Sektor Abfall und Abwasser (CRF 5)

Dieser Sektor besteht vor allem aus Emissionen aus Abfalldeponien, welche der
Zersetzung organischer Materialien in Deponien entstehen. Es wird angenommen,
dass der Abfall aus Produktionsprozessen gegenüber den Abfällen aus Konsum
vernachlässigbar sind, weswegen eine Disaggregation auf Grundlage der
Bevölkerung von ABW vorgenommen wird.

EnbG: Bevölkerung

Quelle:

- [Bevölkerung nach Geschlecht in den Gemeinden](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=12411-0001&bypass=true&levelindex=0&levelid=1691507280245#abreadcrumb)

**Dataset: `raw/emissions`**

??? metadata "Metadata"
    ```json
    {
        "Daten Sachsen-Anhalt": "https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf",
        "Datens\u00e4tze Desaggregation": {
            "Industrie": ""
        }
    }
    ```

------------------------------
## BMWK Langfristszenarien

Langfristszenarien des Bundesministerium für Wirtschaft und Klimaschutz, Daten
auf Deutschlandebene.

Die Daten wurden über den
[Szenario Explorer](https://langfristszenarien.de/enertile-explorer-de/szenario-explorer/)
abgerufen.

### Verwendete Szenarien

- **T45-Strom:** Stromfokussiertes Szenario aus den T45-Szenarien aus 2023, die
  Wege zur Treibhausgasneutralität bis 2045 unter Einhaltung aktueller
  politischer Vorgaben erreichen. Die Daten dieses Szenarios werden als
  Grundlage für das Zielszenario in der Region verwendet.
- **TN-Strom:** Stromfokussiertes Szenario aus den TN-Szenarien aus 2021, die
  unterschiedliche Pfade für Deutschland mit dem Ziel treibhausgasneutral bis
  2050 zu werden. Die Daten dieses Szenarios werden als Grundlage für den
  Status quo verwendet (Ausnahme: Erzeugung Wärmenetze, hier wurden manuell
  Daten für 2021 ergänzt).

### Daten

#### T45-Strom

| Datensatz                                      | Quelle                                                                                                                                                                                                                                                               | Datei                                                     |
|------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| Gebäude: Haushalte und GHD Energiebedarf       | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/51944/21559a9532131c061668bf0751e519e3)                                                                                                                                                            | `T45-Strom_buildings_heating_demand_by_carrier.csv`       |
| Gebäude: Anzahl der Heizungen nach Technologie | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/51944/21559a9532131c061668bf0751e519e3)                                                                                                                                                            | `T45-Strom_buildings_heating_structure_by_technology.csv` |
| GHD Energieträger                              | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798)                                                                                                                                                            | `T45-Strom_cts_demand.csv`                                |
| Haushalte Energieträger                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798)                                                                                                                                                            | `T45-Strom_hh_demand.csv`                                 |
| Industrie Energiebedarf                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52612/9de48084ac2d54c418daaf02a6ee26e0)                                                                                                                                                            | `T45-Strom_ind_demand.csv`                                |
| Stromsystem Deutschland Leistung               | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/48766/5c11999a03c547e04e73d61e4b5fc633)                                                                                                                                                            | `T45-Strom_electricity_installed_power.csv`               |
| Erzeugung Wärmenetze Deutschland               | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/49949/cf898070daec6a4e613dc889927a5feb), [Link2](https://static.agora-energiewende.de/fileadmin/Projekte/2022/2022-11_DE_Large_Scale_Heatpumps/A-EW_293_Rollout_Grosswaermepumpen_WEB.pdf) (S. 37) | `T45-Strom_Generation_Heatgrids_Germany.csv`              |

#### TN-Strom

| Datensatz                                      | Quelle                                                                                                    | Datei                                                    |
|------------------------------------------------|-----------------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| Gebäude: Haushalte und GHD Energiebedarf       | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8198/698cee83d667a2f44fdea7e78ee799a2)  | `TN-Strom_buildings_heating_demand_by_carrier.csv`       |
| Gebäude: Anzahl der Heizungen nach Technologie | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8198/698cee83d667a2f44fdea7e78ee799a2)  | `TN-Strom_buildings_heating_structure_by_technology.csv` |
| GHD Energieträger                              | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8660/ae5a14ff0c320cbd31c5eeff2ede54ba)  | `TN-Strom_cts_demand.csv`                                |
| Haushalte Energieträger                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8660/ae5a14ff0c320cbd31c5eeff2ede54ba)  | `TN-Strom_hh_demand.csv`                                 |
| Industrie Energiebedarf                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/29085/084bd7f45f40d31fd53341e6a94f532c) | `TN-Strom_ind_demand.csv`                                |

**Dataset: `raw/bmwk_long_term_scenarios`**

??? metadata "Metadata"
    ```json
    {
        "name": "bmwk_long_term_scenarios",
        "title": "BMWK Langfristszenarien",
        "id": "bmwk_long_term_scenarios",
        "description": "Langfristszenarien des Bundesministerium f\u00fcr Wirtschaft und Klimaschutz, Daten auf Deutschlandebene.",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "BMWK",
            "Langfristszenario",
            "T45-Strom",
            "TN-Strom"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": null,
            "timeseries": [
                {
                    "start": null,
                    "end": null,
                    "resolution": null,
                    "alignment": null,
                    "aggregationType": null
                },
                {
                    "start": null,
                    "end": null,
                    "resolution": null,
                    "alignment": null,
                    "aggregationType": null
                }
            ]
        },
        "sources": [
            {
                "title": "BMWK Langfristszenarien",
                "description": "Langfristszenarien des Bundesministerium f\u00fcr Wirtschaft und Klimaschutz, Daten auf Deutschlandebene.",
                "path": "https://langfristszenarien.de/enertile-explorer-de/szenario-explorer/",
                "licenses": null
            },
            {
                "title": null,
                "description": null,
                "path": null,
                "licenses": [
                    {
                        "name": null,
                        "title": null,
                        "path": null,
                        "instruction": null,
                        "attribution": null
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": null,
                "title": null,
                "path": null,
                "instruction": null,
                "attribution": null
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-15",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.6.0",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## AGEB – Anwendungsbilanzen für die Endenergiesektoren 2011 bis 2021

Detaillierte Anwendungsbilanzen der Endenergiesektoren für 2020 und 2021 sowie
zusammenfassende Zeitreihen zum Endenergieverbrauch nach Energieträgern und
Anwendungszwecken für Jahre von 2011 bis 2021 der AG Energiebilanzen.

**Dataset: `raw/ageb_energy_balance`**

??? metadata "Metadata"
    ```json
    {
        "name": "ageb_energy_balance",
        "title": "AGEB \u2013 Anwendungsbilanzen f\u00fcr die Endenergiesektoren 2011 bis 2021",
        "id": "ageb_energy_balance",
        "description": "Detaillierte Anwendungsbilanzen der Endenergiesektoren f\u00fcr 2020 und 2021 sowie zusammenfassende Zeitreihen zum Endenergieverbrauch nach Energietr\u00e4gern und Anwendungszwecken f\u00fcr Jahre von 2011 bis 2021",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Endenergiesektoren",
            "Anwendungsbilanzen",
            "energy-balance"
        ],
        "publicationDate": "2022-12-01",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-12-01",
            "timeseries": null
        },
        "sources": [
            {
                "title": "AGEB \u2013 Anwendungsbilanzen f\u00fcr die Endenergiesektoren 2011 bis 2021",
                "description": "Detaillierte Anwendungsbilanzen der Endenergiesektoren f\u00fcr 2020 und 2021 sowie zusammenfassende Zeitreihen zum Endenergieverbrauch nach Energietr\u00e4gern und Anwendungszwecken f\u00fcr Jahre von 2011 bis 2021",
                "path": "https://ag-energiebilanzen.de/daten-und-fakten/anwendungsbilanzen/",
                "licenses": null
            }
        ],
        "licenses": null,
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-15",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": null,
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## sEEnergies Pan-European Thermal Atlas 5.2 (Peta5)

Wärmebedarf für Europa 2015 in GJ (1ha Auflösung) für

- Haushalte: Raumwärme und Warmwasser
- GHD: Raumwärme, Warmwasser und Prozesswärme

Die Daten können auf der
[Projektseite](https://s-eenergies-open-data-euf.hub.arcgis.com)
eingesehen werden.

### Haushalte

Abgerufen mittels

```commandline
wget -O Peta5_0_1_HD_res.zip https://arcgis.com/sharing/rest/content/items/d7d18b63250240a49eb81db972aa573e/data
```

### GHD und Industrie

Abgerufen mittels

```commandline
wget -O Peta5_0_1_HD_ser.zip https://arcgis.com/sharing/rest/content/items/52ff5e02111142459ed5c2fe3d80b3a0/data
```

**Dataset: `raw/seenergies_peta5`**

??? metadata "Metadata"
    ```json
    {
        "name": "seenergies_peta5",
        "title": "sEEnergies Pan-European Thermal Atlas 5.2 (Peta5)",
        "id": "seenergies_peta5",
        "description": "W\u00e4rmebedarf f\u00fcr Europa 2015 in GJ (1ha Aufl\u00f6sung)",
        "language": [
            "en-GB"
        ],
        "subject": null,
        "keywords": [
            "European",
            "Photovoltaic",
            "Fl\u00e4chenpotential",
            "Aufdachanlagen"
        ],
        "publicationDate": "2022-01-01",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Europe",
            "extent": "Europe",
            "resolution": "1 ha"
        },
        "temporal": {
            "referenceDate": "2022-01-01",
            "timeseries": null
        },
        "sources": [
            {
                "title": "sEEnergies Pan-European Thermal Atlas 5.2 (Peta5)",
                "description": "W\u00e4rmebedarf f\u00fcr Europa 2015 in GJ (1ha Aufl\u00f6sung)",
                "path": "https://www.seenergies.eu/peta5/",
                "licenses": [
                    {
                        "name": null,
                        "title": null,
                        "path": null,
                        "instruction": "The data provided is indiative and for research purpose only",
                        "attribution": "\u00a9 Flensburg, Halmstad and Aalborg Universities 2022"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": null,
                "title": null,
                "path": null,
                "instruction": "The data provided is indiative and for research purpose only",
                "attribution": "\u00a9 Flensburg, Halmstad and Aalborg Universities 2022"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Sozialversicherungspflichtig Beschäftigte und Betriebe

Gemeindedaten der sozialversicherungspflichtig Beschäftigten am 30.06.2022 nach
Wohn- und Arbeitsort - Deutschland, Länder, Kreise und Gemeinden (Jahreszahlen)
der Bundesagentur für Arbeit.

**Dataset: `raw/ba_employment`**

??? metadata "Metadata"
    ```json
    {
        "name": "ba_employment",
        "title": "Gemeindedaten der sozialversicherungspflichtig Besch\u00e4ftigten nach Wohn- und Arbeitsort",
        "id": "ba_employment",
        "description": "Zahl der soziaversicherungspflichtig Besch\u00e4ftigten nach: Wohnort, Personengruppen, Arbeitsort, Wohnort gleich Arbeitsort, Einpendeler, Auspendler, Zahl der Betriebe",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Gemeindedaten",
            "sozialversicherungspflichtig",
            "Besch\u00e4ftigte",
            "Wohnort",
            "Arbeitsort"
        ],
        "publicationDate": "2023-01-16",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-06-23",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Gemeindedaten der sozialversicherungspflichtig Besch\u00e4ftigten nach Wohn- und Arbeitsort",
                "description": "Zahl der soziaversicherungspflichtig Besch\u00e4ftigten nach: Wohnort, Personengruppen, Arbeitsort, Wohnort gleich Arbeitsort, Einpendeler, Auspendler, Zahl der Betriebe",
                "path": "https://statistik.arbeitsagentur.de/SiteGlobals/Forms/Suche/Einzelheftsuche_Formular.html?nn=15024&topic_f=beschaeftigung-sozbe-gemband",
                "licenses": [
                    {
                        "name": null,
                        "title": null,
                        "path": "https://statistik.arbeitsagentur.de/DE/Statischer-Content/Servicebereich-Navigation/Bezugsbedingungen.html?nn=6654",
                        "instruction": "Sie k\u00f6nnen Informationen speichern, (auch auszugsweise) mit Quellenangabe weitergeben, vervielf\u00e4ltigen und verbreiten. Die Inhalte d\u00fcrfen nicht ver\u00e4ndert oder verf\u00e4lscht werden. Eigene Berechnungen sind erlaubt, jedoch als solche kenntlich zu machen. Im Falle einer Zug\u00e4nglichmachung im Internet soll dies in Form einer Verlinkung auf die Homepage der Statistik der Bundesagentur f\u00fcr Arbeit erfolgen. Die Nutzung der Inhalte f\u00fcr gewerbliche Zwecke, ausgenommen Presse, Rundfunk und Fernsehen und wissenschaftliche Publikationen, bedarf der Genehmigung durch die Statistik der Bundesagentur f\u00fcr Arbeit.",
                        "attribution": "\u00a9 Statistik der Bundesagentur f\u00fcr Arbeit"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": null,
                "title": null,
                "path": "https://statistik.arbeitsagentur.de/DE/Statischer-Content/Servicebereich-Navigation/Bezugsbedingungen.html?nn=6654",
                "instruction": "Sie k\u00f6nnen Informationen speichern, (auch auszugsweise) mit Quellenangabe weitergeben, vervielf\u00e4ltigen und verbreiten. Die Inhalte d\u00fcrfen nicht ver\u00e4ndert oder verf\u00e4lscht werden. Eigene Berechnungen sind erlaubt, jedoch als solche kenntlich zu machen. Im Falle einer Zug\u00e4nglichmachung im Internet soll dies in Form einer Verlinkung auf die Homepage der Statistik der Bundesagentur f\u00fcr Arbeit erfolgen. Die Nutzung der Inhalte f\u00fcr gewerbliche Zwecke, ausgenommen Presse, Rundfunk und Fernsehen und wissenschaftliche Publikationen, bedarf der Genehmigung durch die Statistik der Bundesagentur f\u00fcr Arbeit.",
                "attribution": "\u00a9 Statistik der Bundesagentur f\u00fcr Arbeit"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-15",
                "object": "metadata",
                "comment": "Create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Bevölkerung

Einwohnerzahl nach Gemeinden des Statistischen Bundesamts.

**Dataset: `raw/destatis_gv`**

??? metadata "Metadata"
    ```json
    {
        "name": "destatis_gv",
        "title": "Adminstratives Gemeinndeverzeichnis",
        "id": "destatis_gv",
        "description": "Alle politisch selbst\u00e4ndigen Gemeinden mit ausgew\u00e4hlten Merkmalen am 31.12.2022 ",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "destatis",
            "gemeindeverzeichnis"
        ],
        "publicationDate": "2023-01-12",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": "2022-02-14",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Statistisches Bundesamt",
                "description": "Alle politisch selbst\u00e4ndigen Gemeineden mit ausgew\u00e4hlten Merkmalen am 31.12.2022 (4.Quartal)",
                "path": "https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/Archiv/GVAuszugQ/AuszugGV4QAktuell.html",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Statistisches Bundesamt (Destatis), 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Statistisches Bundesamt (Destatis), 2023"
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Regionalstatistik (GENESIS)

Enthält folgende Datensätze der statistischen Ämter des Bundes und der Länder:

### Energieverwendung der Betriebe im Verarbeitenden Gewerbe (43531-01-02-4)

Jahreserhebung ü. die Energieverwendung der Betriebe im verarbeitendem Gewerbe.

Der Datensatz umfasst:

- Betriebe des Verarbeitenden Gewerbes sowie des Bergbaus
und der Gewinnung von Steinen und Erden von Unternehmen des
Produzierenden Gewerbes mit im Allgemeinen 20 und mehr
Beschäftigten.
- Betriebe des Verarbeitenden Gewerbes sowie des Bergbaus
und  der Gewinnung von Steinen und Erden mit im Allgemeinen
20 und mehr Beschäftigten von Unternehmen der übrigen
Wirtschaftsbereiche.
Die Berichterstattung schließt Verarbeitende Betriebe des
Handwerks ein.
Bei 7 Wirtschaftszweigen gilt eine Abschneidegrenze von 10
Beschäftigten. Die Merkmalswerte beziehen sich auf den
gesamten Betrieb, schließen damit die nicht produzierenden
Betriebsteile mit ein.
Maßgebend für die Zuordnung ist ab 2008 die „Klassifikation
der Wirtschaftszweige, Ausgabe 2008 (WZ 2008)“, und zwar
die Abschnitte B und C.

- Datei: `43531-01-02-4.xlsx`
- Stand: 2021

### Betriebe, tätige Personen, Bruttoentgelte (42111-01-04-5)

Jahreserhebung ü. Betriebe, tätige Personen und Bruttoentgelte der Betriebe im
verarbeitendem Gewerbe.

Der Datensatz umfasst:

- Sämtliche Betriebe des Wirtschaftsbereiches Verarbeitendes
Gewerbe sowie Bergbau und Gewinnung von Steinen und Erden,
wenn diese Betriebe zu Unternehmen des Bereiches
Verarbeitendes Gewerbe sowie Bergbau und Gewinnung von
Steinen und Erden gehören und in diesen Unternehmen
mindestens 20 Personen tätig sind;
- die Betriebe des Wirtschaftsbereiches Verarbeitendes
Gewerbe sowie Bergbau und Gewinnung von Steinen und Erden
mit mindestens 20 tätigen Personen, sofern diese Betriebe
zu Unternehmen gehören, deren wirtschaftlicher Schwerpunkt
außerhalb des Bereiches Verarbeitendes Gewerbe sowie
Bergbau und Gewinnung von Steinen und Erden liegt.
Bei 7 kleinbetrieblich strukturierten Branchen gilt eine
untere Erfassungsgrenze von 10 tätigen Personen.
Die Auswahl erfolgt jeweils nach dem Beschäftigtenstand Ende
September des Vorjahres. Die ausgewiesene Beschäftigtenzahl
betrifft dagegen die von Ende September des Berichtsjahres.
Die Merkmalswerte beziehen sich auf den gesamten Betrieb,
schließen damit die nicht produzierenden Betriebsteile mit
ein.
Maßgebend für die Zuordnung ist ab 2009 die „Klassifikation
der Wirtschaftszweige, Ausgabe 2008 (WZ 2008)“, und zwar
die Abschnitte B und C.

- Datei: `42111-01-04-5.xlsx`
- Stand: 30.09.2021

### Gebäude mit Wohnraum nach Heizungsart (31211-04-01-5-B)

Zensus 2011: Gebäude mit Wohnraum nach Heizungsart

- Datei: `31211-04-01-5-B.xlsx`
- Stand: 09.05.2011

### Gebäude mit Wohnraum nach Heizungsart (31231-02-01-5)

Bestand an Wohngebäuden und Wohnungen in Wohn- und Nichtwohngebäuden -
Fortschreibung auf Basis der endgültigen Ergebnisse der Gebäude- und
Wohnungszählung 2011 (Zensus 2011).

- Datei: `31231-02-01-5.xlsx`
- Stand: 31.12.2021

**Dataset: `raw/regiostat`**

??? metadata "Metadata"
    ```json
    {
        "name": "regiostat",
        "title": "Regionalstatistik (GENESIS)",
        "id": "regiostat",
        "description": "Energieverwendung der Betriebe im Verarbeitenden Gewerbe (43531-01-02-4), Betriebe, t\u00e4tige Personen, Bruttoentgelte (42111-01-04-5), Geb\u00e4ude mit Wohnraum nach Heizungsart (31211-04-01-5-B), Geb\u00e4ude mit Wohnraum nach Heizungsart (31231-02-01-5)",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Regionalstatistik",
            "Energieverwendung",
            "verarbeitendes Gewerbe",
            "t\u00e4tige Personen",
            "Bruttoentgelte",
            "Geb\u00e4ude",
            "Heizungsart"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2021-01-01",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Regionalstatistik (GENESIS)",
                "description": "Energieverwendung der Betriebe im Verarbeitenden Gewerbe (43531-01-02-4), Betriebe, t\u00e4tige Personen, Bruttoentgelte (42111-01-04-5), Geb\u00e4ude mit Wohnraum nach Heizungsart (31211-04-01-5-B), Geb\u00e4ude mit Wohnraum nach Heizungsart (31231-02-01-5)",
                "path": [
                    "https://www.regionalstatistik.de/genesis//online?operation=table&code=43531-01-02-4",
                    "https://www.regionalstatistik.de/genesis//online?operation=table&code=42111-01-04-5"
                ],
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9  Statistische \u00c4mter des Bundes und der L\u00e4nder, 2023"
                    }
                ]
            },
            {
                "title": null,
                "description": null,
                "path": null,
                "licenses": [
                    {
                        "name": null,
                        "title": null,
                        "path": null,
                        "instruction": null,
                        "attribution": null
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": null,
                "title": null,
                "path": null,
                "instruction": null,
                "attribution": null
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-15",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Geodaten PV- und Windflächenrechner

Geodaten aus dem [PV- und Windflächenrechner](https://www.agora-energiewende.de/service/pv-und-windflaechenrechner/).

Mehr Informationen:

- [Begleitdokument](https://zenodo.org/record/6794558)
- [Geodaten Potenzialflächen](https://zenodo.org/record/6728382)

Enthält:

- Geodaten
- Metadaten
- App-Datapackage

**Dataset: `raw/rli_pv_wfr`**

??? metadata "Metadata"
    ```json
    {
        "name": "rli_pv_wfr",
        "title": "Geodaten PV- und Windfl\u00e4chenrechner",
        "id": "rli_pv_wfr",
        "description": "Geodaten aus dem PV- und Windfl\u00e4chenrechner",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Geodaten",
            "PV-Fl\u00e4chenrechner",
            "Windfl\u00e4chenrechner",
            "Potentialfl\u00e4chen"
        ],
        "publicationDate": "2022-06-05",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-06-05",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Geodaten PV- und Windfl\u00e4chenrechner",
                "description": "Geodaten aus dem PV- und Windfl\u00e4chenrechner",
                "path": "https://zenodo.org/record/6728382",
                "licenses": null
            },
            {
                "title": null,
                "description": null,
                "path": null,
                "licenses": [
                    {
                        "name": "CC BY-NC 4.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://creativecommons.org/licenses/by-nc/4.0/",
                        "instruction": "you are free to copy, redistribute and adapt them for non-commercial purposes, provided you give appropriate credit. Note that the data is made available as-is and without warranty. We cannot guarantee its accuracy, and accept no responsibility for any liability arising from its use. You are advised to examine the quality of the data for your intended purposes, and to consult the publications linked on this page.",
                        "attribution": "\u00a9 Reiner Lemoine INstitut, 2022"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "CC BY-NC 4.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://creativecommons.org/licenses/by-nc/4.0/",
                "instruction": "you are free to copy, redistribute and adapt them for non-commercial purposes, provided you give appropriate credit. Note that the data is made available as-is and without warranty. We cannot guarantee its accuracy, and accept no responsibility for any liability arising from its use. You are advised to examine the quality of the data for your intended purposes, and to consult the publications linked on this page.",
                "attribution": "\u00a9 Reiner Lemoine INstitut, 2022"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-15",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Regionalplan Anhalt-Bitterfeld-Wittenberg

Geodatensätze aus Teilplänen Wind 2018 und 2027 der Regionalen
Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg.

### Sachlicher Teilplan Wind 2018

Geodaten aus rechtskräftigem
[Sachlichen Teilplan Wind 2018](https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2018/).

> Im Sachlichen Teilplan "Nutzung der Windenergie in der Planungsregion
> Anhalt-Bitterfeld-Wittenberg" vom 30.05.2018 werden 22 Vorranggebiete für die
> Nutzung der Windenergie mit der Wirkung von Eignungsgebieten festgelegt. Sie
> dienen der raumordnerischen Steuerung der Errichtung von raumbedeutsamen
> Windenergieanlagen in Konzentrationszonen.
>
> Die oberste Landesentwicklungsbehörde hat am 01.08.2018 die Genehmigung
> erteilt. Mit Bekanntmachung der Genehmigung tritt der Sachliche Teilplan in
> Kraft.

Dateien:

- Vorrang-/Eignungsgebiete: `stp_2018_vreg.gpkg`
  ([Quelle](https://gis.planungsregion-abw.de/geoserver/stp_wind2018/ows?SERVICE=WFS&REQUEST=GetCapabilities))

### Sachlicher Teilplan Wind 2027

Geodaten aus Planentwurf des
[Sachlichen Teilplan Wind 2027](https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2027/).

> Die Regionalversammlung hat am 03.03.2023 beschlossen, den Sachlichen
> Teilplan "Windenergie 2027 in der Planungsregion Anhalt-Bitterfeld-Wittenberg"
> aufzustellen und mit der Bekanntgabe der Allgemeinen Planungsabsicht die
> beabsichtigten Auswahlkriterien und mögliche Gebietskulisse der Vorranggebiete
> für die Nutzung der Windenergie bzw. für Repowering von Windenergieanlagen
> vorzustellen.

Dateien:

- Suchräume: `stp_2027_suchraum.gpkg` (Quelle: RPG ABW)
- Planabsicht Vorranggebiete: `stp_2027_ideen_vr.gpkg` (Quelle: RPG ABW)
- Planabsicht Repoweringgebiete: `stp_2027_ideen_repower.gpkg` (Quelle: RPG ABW)

**Dataset: `raw/rpg_abw_regional_plan`**

??? metadata "Metadata"
    ```json
    {
        "name": "rpg_abw_regional_plan",
        "title": "Regionalplan Anhalt-Bitterfeld-Wittenberg",
        "id": "rpg_abw_regional_plan",
        "description": "Geodatens\u00e4tze aus Teilpl\u00e4nen Wind 2018 und 2027 der Regionalen Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg.",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Geodatens\u00e4tze",
            "Teilpl\u00e4ne",
            "PLanungsgemeinschaft"
        ],
        "publicationDate": "2018-05-30",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "ABW",
            "extent": "ABW",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2023-03-03",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Regionalplan Anhalt-Bitterfeld-Wittenberg",
                "description": "Geodatens\u00e4tze aus Teilpl\u00e4nen Wind 2018 und 2027 der Regionalen Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg.",
                "path": [
                    "https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2018/",
                    "https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2027/"
                ],
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Regionale Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg / Jahr"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Regionale Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg / Jahr"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Marktstammdatenregister Datenkorrektur PV

Überprüfung und manuelle Datenkorrektur der Photovoltaikanlagen aus dem
prozessierten Marktstammdatenregister (Datensatz:
[bnetza_mastr](../../digipipe/store/raw/bnetza_mastr/dataset.md)).

### Plausibiltätsprüfung

Um grobe Fehler herauszufiltern wird überprüft, ob

- Anlage in Betrieb ist (status = "In Betrieb"),
- Anlage Strom produziert,
- Brutto- und Nettokapazität plausibel sind und
- die Kategorisierung, d.h. Zuordnung eine PV-Anlage zu Freifläche oder Dach,
  plausibel ist (manuelle, visuelle Prüfung von geolokalisierten
  PV-Aufdachanlagen anhand von
  [Orthofotos](https://www.geodatenportal.sachsen-anhalt.de/wss/service/ST_LVermGeo_DOP_WMS_OpenData/guest))

### Dateien

- Korrektur Freiflächenanlagen `bnetza_mastr_pv_ground_region_correction.ods`
- Korrektur Aufdachanlagen `bnetza_mastr_pv_roof_region_correction.ods`

mit Spalten:

- _mastr_id_: ID aus dem MaStR
- _reason_: Fehler (wrong_type, wrong_position)
- _wrong_attr_: Fehlerhaftes Attribut
- _correction_: Korrigierter Attributwert (None, wenn Korrektur nicht möglich).
  Korrigierte Geometrien liegen in EPSG:3035 vor.

**Dataset: `raw/bnetza_mastr_correction_region`**

??? metadata "Metadata"
    ```json
    {
        "name": "bnetza_mastr_correction",
        "title": "Marktstammdatenregisterdaten - Manuelle Korrektur",
        "id": "bnetza_mastr",
        "description": "Daten aus dem Marktstammdatenregister der Bundesnetzagentur",
        "language": [
            "en-GB",
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Markstammdatenregister",
            "openmastr",
            "mastr"
        ],
        "publicationDate": "2022-12-19",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-12-19",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Marktstammdatenregister",
                "description": "Marktstammdatenregister der Bundesnetzagentur Deutschland",
                "path": "https://www.marktstammdatenregister.de/MaStR/Datendownload",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0",
                        "path": "http://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Marktstammdatenregister 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0?",
                "path": "http://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Marktstammdatenregister 2023"
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "Create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Erzeugungsanlagen aus Marktstammdatenregister

Ereugungsanlagen aus dem Markstammdatenregister, das mit dem Tool
[open-mastr](https://github.com/OpenEnergyPlatform/open-MaStR) erstellt und
abgelegt wurde. Die Daten wurden folgendermaßen erstellt:
```
from open_mastr import Mastr
db = Mastr()
db.download("bulk")
db.to_csv(None)  # (None for all data)
```

Die abgelegten CSV-Dateien (alle Tabellen) wurden um einen benutzerdefinierten
Export von Speichereinheiten mit
`sqlite3 -header -csv -separator "," open-mastr.db "select * from storage_units;" > bnetza_mastr_storage_unit_raw.csv`
erweitert. Anschließend wurden alle Dateien komprimiert.

Das Marktstammdatenregister (MaStR) ist ein deutsches Register, welches von der
Bundesnetzagentur (BNetza) bereitgestellt wird und alle in Deutschland
befindlichen Strom- und Gasanlagen erfasst.

**Dataset: `raw/bnetza_mastr`**

??? metadata "Metadata"
    ```json
    {
        "name": "bnetza_mastr",
        "title": "Marktstammdatenregisterdaten",
        "id": "bnetza_mastr",
        "description": "Daten aus dem Marktstammdatenregister der Bundesnetzagentur",
        "language": [
            "en-GB",
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Markstammdatenregister",
            "openmastr",
            "mastr"
        ],
        "publicationDate": "2022-12-19",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-12-19",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Marktstammdatenregister",
                "description": "Marktstammdatenregister der Bundesnetzagentur Deutschland",
                "path": "https://www.marktstammdatenregister.de/MaStR/Datendownload",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0",
                        "path": "http://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Marktstammdatenregister 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0?",
                "path": "http://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Marktstammdatenregister 2023"
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "Create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Energiedaten Sachsen-Anhalt

Datensätze zur Energie- und Wasserversorgung des Statistischen Landesamtes
Sachsen-Anhalt.

### Daten

Stromverbrauch der Industriebetriebe nach Kreisen 2003-2021 in MWh

- [Quelle](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c206986)

**Dataset: `raw/stala_st_energy`**

??? metadata "Metadata"
    ```json
    {
        "name": "stala_st_energy",
        "title": "Energiedaten Sachsen-Anhalt",
        "id": "stala_st_energy",
        "description": "Datens\u00e4tze zur Energie- und Wasserversorgung des Statistischen Landesamtes Sachsen-Anhalt.",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Energiedaten",
            "Energieversorgung",
            "Wasserversorgung"
        ],
        "publicationDate": "2022-01-01",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Sachsen-Anhalt",
            "extent": "Sachsen-Anhalt",
            "resolution": "NUTS-3"
        },
        "temporal": {
            "referenceDate": "2022-01-01",
            "timeseries": [
                {
                    "start": "2003-01-01T00:00+01",
                    "end": "2021-12-31T23:00+01",
                    "resolution": "1 a",
                    "alignment": "left"
                }
            ]
        },
        "sources": [
            {
                "title": "Energiedaten Sachsen-Anhalt",
                "description": "Datens\u00e4tze zur Energie- und Wasserversorgung des Statistischen Landesamtes Sachsen-Anhalt.",
                "path": "https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c206986",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Statistisches Landesamt Sachsen-Anhalt, Halle (Saale)."
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Statistisches Landesamt Sachsen-Anhalt, Halle (Saale)."
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## DemandRegio

Regionalisierte Bevölkerungsprognose, Haushalte sowie Strom- und Gasbedarfe
inkl. Zeitreihen auf Landkreisebene.

Die Daten wurden abgerufen mit einer
[modifizierten Version des DemandRegio disaggregators](https://github.com/nesnoj/disaggregator),
in der softwareseitige, jedoch keine methodischen Änderungen vorgenommen wurden.

Der disaggregator basiert auf Daten bis 2017, anschließende Jahre werden
fortgeschrieben.

Weitere Informationen zum Projekt DemandRegio:

- [Abschlussbericht](https://www.ffe.de/wp-content/uploads/2020/10/DemandRegio_Abschlussbericht.pdf)
- [Abschlussworkshop](https://www.tu.berlin/er/forschung/projekte/demandregio-2)

Die erzeugten Rohdaten wie unten beschrieben wurden mittels
[API](http://opendata.ffe.de:4000/) abgerufen. Diese können alternativ direkt
vom [OpenData-Portal der FfE](https://opendata.ffe.de/project/demandregio/)
bezogen werden.

Verwendetes Wetterjahr für Gasbedarfszeitreihen: 2011

**Installation (in separater venv):**

```commandline
pip install disaggregator@git+https://github.com/nesnoj/disaggregator.git#egg=disaggregator
```

### Details zum Datenabruf

#### Bevölkerung

Bevölkerung (Summe) und Bevölkerung je Haushaltsgröße (1, 2, 3, 4, 5, >5) je
NUTS3.

Jahre:

- Bevölkerung bis 2017 historische Werte
- Bevölkerung ab 2018 prognostizierte Werte basierend auf der 14. koordinierten
  Bevölkerungsvorausberechnung der Statistischen Ämter von Bund und Ländern.
- Haushalte nur 2011

```python
import pandas as pd
from disaggregator import data

## Population
dr_hh_population = pd.DataFrame()
for year in [2010, 2015, 2017, 2020, 2021, 2022, 2025, 2030, 2035, 2040, 2045]:
    dr_hh_population[year] = round(data.population(year=year)).astype(int)

dr_hh_population.to_csv("dr_hh_population.csv")

## Households
data.households_per_size().to_csv("dr_hh_households_2011.csv")
```

#### Haushalte: Strom

Bedarfe und SLP-Zeitreihen je NUTS3 mit Bottom-Up-Methode nach Haushaltsgröße.

Jahre:

- 2017: Letzte verfügbare Daten
- 2022: Status quo, Fortschreibung mit Berücksichtigung Demografie und
  Wanderung
- 2035: Fortschreibungsjahr mit Berücksichtigung Demografie und Wanderung
- 2045: Fortschreibungsjahr

```python
from disaggregator import spatial, temporal

## Consumption
spatial.disagg_households_power(
    by="households",
    weight_by_income=True,
    year=2022,
    scale_by_pop=True,
).to_csv(f"dr_hh_power_demand_2022.csv")

## Timeseries
temporal.disagg_temporal_power_housholds_slp(
    use_nuts3code=True,
    by="households",
    weight_by_income=True,
    year=2022,
    scale_by_pop=True,
).to_csv(f"dr_hh_power_timeseries_2022.csv")
```

#### Haushalte: Gas

Zeitreihen je NUTS3

```python
from disaggregator import temporal

## Timeseries
temporal.disagg_temporal_gas_households(
    use_nuts3code=True,
    how='top-down',
    year=2011,
).to_csv(f"dr_hh_gas_timeseries_2011.csv")
```


#### GHD und Industrie: Strom

Bedarfe und Zeitreihen je NUTS3:

- Bedarfe: Je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: Für alle WZ bedarfsgewichtet aggregiert, Einzelprofile basieren
  je nach WZ auf gemessenen oder SLP inkl. Wanderung
- Letzte verfügbare Daten aus 2017, Fortschreibung für 2022 mit
  Berücksichtigung Beschäftigte und Effizienzgewinne

```python
from disaggregator import spatial, temporal

########
## CTS #
########

## Consumption
spatial.disagg_CTS_industry(
    sector='CTS',
    source='power',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_cts_power_demand_2022.csv")
## Timeseries
temporal.disagg_temporal_power_CTS(
    detailed=False,
    use_nuts3code=True,
    year=2022,
).to_csv("dr_cts_power_timeseries_2022.csv")

#############
## Industry #
#############

## Consumption
spatial.disagg_CTS_industry(
    sector='industry',
    source='power',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_ind_power_demand_2022.csv")
## Timeseries
temporal.disagg_temporal_industry(
    source="power",
    detailed=False,
    use_nuts3code=True,
    no_self_gen=False,
    year=2022,
).to_csv("dr_ind_power_timeseries_2022.csv")
```

#### GHD: Gas

Zeitreihen je NUTS3 für alle WZ bedarfsgewichtet aggregiert, Einzelprofile
basieren je nach WZ auf gemessenen oder SLP inkl. Wanderung. Letzte verfügbare
Daten aus 2017, Fortschreibung für 2022 mit Berücksichtigung Beschäftigte und
Effizienzgewinne.

```python
from disaggregator import spatial, temporal

## Timeseries
x=temporal.disagg_temporal_gas_CTS(
    detailed=False,
    use_nuts3code=True,
    year=2011,
).to_csv("dr_cts_gas_timeseries_2011.csv")
```

#### Industrie: Gas

Bedarfe und Zeitreihen je NUTS3:

- Bedarfe: Je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: Für alle WZ bedarfsgewichtet aggregiert, Einzelprofile basieren
  je nach WZ auf gemessenen oder SLP inkl. Wanderung
- Letzte verfügbare Daten aus 2017, Fortschreibung für 2022 mit
  Berücksichtigung Beschäftigte und Effizienzgewinne

```python
from disaggregator import spatial, temporal

## Consumption
spatial.disagg_CTS_industry(
    sector='industry',
    source='gas',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_ind_gas_demand_2022.csv")
## Timeseries
x=temporal.disagg_temporal_industry(
    source="gas",
    detailed=False,
    use_nuts3code=True,
    no_self_gen=False,
    year=2011,
).to_csv("dr_ind_gas_timeseries_2011.csv")
```

**Dataset: `raw/demandregio`**

??? metadata "Metadata"
    ```json
    {
        "name": "demandregio",
        "title": "DemandRegio",
        "id": "demandregio",
        "description": "Regionalisierte Bev\u00f6lkerungsprognose, Haushalte sowie Strom- und Gasbedarfe inkl. Zeitreihen auf Landkreisebene.",
        "language": [
            "en-GB",
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Bev\u00f6lkerung",
            "Bev\u00f6lkerungsprognose",
            "Strombedarf",
            "Gasbedarf",
            "Haushalte",
            "disaggregator",
            "disaggregation"
        ],
        "publicationDate": "2020-09-30",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": "NUTS-3"
        },
        "temporal": {
            "referenceDate": "2022-01-01",
            "timeseries": [
                {
                    "start": "2022-01-01T00:00+01",
                    "end": "2022-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": null,
                    "aggregationType": "sum"
                },
                {
                    "start": "2035-01-01T00:00+01",
                    "end": "2035-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": null,
                    "aggregationType": "sum"
                },
                {
                    "start": "2045-01-01T00:00+01",
                    "end": "2045-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": null,
                    "aggregationType": "sum"
                },
                {
                    "start": "2011-01-01T00:00+01",
                    "end": "2011-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": null,
                    "aggregationType": "sum"
                },
                {
                    "start": "2022-01-01T00:00+01",
                    "end": "2022-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": null,
                    "aggregationType": "sum"
                }
            ]
        },
        "sources": [
            {
                "title": "DemandRegio",
                "description": "Regionalisierte Bev\u00f6lkerungsprognose, Haushalte sowie Strom- und Gasbedarfe inkl. Zeitreihen auf Landkreisebene.",
                "path": "https://github.com/nesnoj/disaggregator/",
                "licenses": [
                    {
                        "name": "CC BY 4.0 DE",
                        "title": "Creative Commons Namensnennung 4.0 Deutschland",
                        "path": "https://creativecommons.org/licenses/by/4.0/deed.de",
                        "instruction": "Sie m\u00fcssen angemessene Urheber- und Rechteangaben machen, einen Link zur Lizenz beif\u00fcgen und angeben, ob \u00c4nderungen vorgenommen wurden. Diese Angaben d\u00fcrfen in jeder angemessenen Art und Weise gemacht werden, allerdings nicht so, dass der Eindruck entsteht, der Lizenzgeber unterst\u00fctze gerade Sie oder Ihre Nutzung besonders.",
                        "attribution": "\u00a9 FZJ, TUB, FfE"
                    },
                    {
                        "name": "GNU FDL",
                        "title": "GNU General Public License v3.0",
                        "path": "https://www.gnu.org/licenses/gpl-3.0.en.html",
                        "instruction": "Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed",
                        "attribution": "\u00a9 FZJ, TUB, FfE"
                    }
                ]
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-08-15",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## OpenStreetMap

OpenStreetMap Datenauszug Deutschland.

Quelle: https://download.geofabrik.de/europe/germany-230101.osm.pbf

Ist nicht Teil des Eingangsdaten-Packages - manueller Download erforderlich.

**Dataset: `raw/osm`**

??? metadata "Metadata"
    ```json
    {
        "name": "openstreetmap",
        "title": "",
        "id": "openstreetmap",
        "description": "OpenStreetMap extract",
        "language": [
            "de-DE",
            "en-GB"
        ],
        "subject": [],
        "keywords": [
            "openstreetmap",
            "osm"
        ],
        "publicationDate": "2023-06-30",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": "2023-06-30",
            "timeseries": []
        },
        "sources": [
            {
                "title": "OpenStreetMap Data Extracts (Geofabrik)",
                "description": "Full data extract of OpenStreetMap data",
                "path": "https://download.geofabrik.de/europe/germany-230101.osm.pbf",
                "licenses": [
                    {
                        "name": "ODbL-1.0",
                        "title": "Open Data Commons Open Database License 1.0",
                        "path": "https://opendatacommons.org/licenses/odbl/1.0/",
                        "instruction": "You are free: To Share, To Create, To Adapt; As long as you: Attribute, Share-Alike, Keep open!",
                        "attribution": "\u00a9 OpenStreetMap contributors"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "ODbL-1.0",
                "title": "Open Data Commons Open Database License 1.0",
                "path": "https://opendatacommons.org/licenses/odbl/1.0/",
                "instruction": "You are free: To Share, To Create, To Adapt; As long as you: Attribute, Share-Alike, Keep open!",
                "attribution": "\u00a9 OpenStreetMap contributors"
            }
        ],
        "contributors": [
            {
                "title": "nesnoj",
                "email": "jonathan.amme@rl-institut.de",
                "date": "2023-06-30",
                "object": "metadata",
                "comment": "Create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "oemetadata_v1.5.1",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## EE-Einspeisezeitreihen

Einspeisezeitreihen für Erneuerbare Energien, normiert auf 1 MW bzw. 1 p.u.
Als Wetterjahr wird 2011 verwendet, siehe
[Szenarien](../../digipipe/store/../../docs/sections/scenarios.md).

### Windenergie

Stündlich aufgelöste Zeitreihe der Windenergie Einspeisung über 1 Jahr auf Basis
von [MaStR](../../digipipe/store/raw/bnetza_mastr/dataset.md) und
[renewables.ninja](http://renewables.ninja).
Auf einen Auflösung auf Gemeindeebene wird verzichtet, da die Differenz der
Produktion der Gemeinden nach renewables.ninja <5 % beträgt.

#### Windenergieanlage (2022)

Für renewables.ninja sind Position (lat, lon), Nennleistung (capacity),
Nabenhöhe und Turbinentyp erforderlich.

##### Position

Hierfür wird aus den Zentroiden der Gemeinden ein räumlicher Mittelwert
anhand des Datensatzes
[bkg_vg250_muns_region](../../digipipe/store/datasets/bkg_vg250_muns_region/dataset.md)
(`bkg_vg250_muns_region.gpkg`) gebildet:

```
import geopandas as gpd
import os.path

def get_position(gdf):
    df = gpd.read_file(gdf)
    points_of_muns = df["geometry"].centroid
    points_of_muns_crs = points_of_muns.to_crs(4326)
    point_df = [
        points_of_muns_crs.y.sum()/len(points_of_muns),
        points_of_muns_crs.x.sum()/len(points_of_muns)
    ]
    return point_df

data_folder = os.path.join("your_data_folder")
muns_gpkg = os.path.join(data_folder, "bkg_vg250_muns_region.gpkg")
center_position = get_position(muns_gpkg)
```

##### Nennleistung

Wird auf 1 MW gesetzt/normiert.

##### Nabenhöhe

Aus dem Datensatz
[bnetza_mastr_wind_region](../../digipipe/store/datasets/bnetza_mastr_wind_region/dataset.md)
(`bnetza_mastr_wind_agg_abw.gpkg`) wird ein Mittelwer von 100 m abgeleitet.

```
import geopandas as gpd

df = gpd.read_file("bnetza_mastr_wind_agg_abw.gpkg")
height = df[["hub_height"]].mean()
```

##### Turbinentyp

Annahme: Innerhalb eines Herstellers sind Leistungskurven sehr ähnlich.
Daher werden zwei größten Hersteller mit jeweiligen häufigsten Turbinentyp
ausgewählt - diese sind Enercon und Vestas mit ca. 70 % und ca. 30%.

```
import geopandas as gpd

df = gpd.read_file("bnetza_mastr_wind_agg_abw.gpkg")
manufacturers = df[
    ["manufacturer_name", "status"]
].groupby("manufacturer_name").count().sort_values(
    by="status", ascending=False
)
```

Häufigste Turbinentypen sind *Enercon E-70* und *Vestas V80*. Daher werden
*Enercon E70 2000* und *Vestas V80 2000* in renewables.ninja ausgewählt.

```
man_1 = manufacturers.index[0]
man_2 = manufacturers.index[1]

type_1 = df[
    ["manufacturer_name", "type_name", "status"]
].where(df["manufacturer_name"] == man_1).groupby(
    "type_name").count().sort_values(by="status", ascending=False)

type_2 = df[
    ["manufacturer_name", "type_name", "status"]
].where(df["manufacturer_name"] == man_2).groupby(
    "type_name").count().sort_values(by="status", ascending=False)
```

#### Raw Data von [renewables.ninja](http://renewables.ninja) API

Es werden zwei Zeitreihen für oben beschriebenen Vergleichsanlagen berechnet:

```
import json
import requests
import pandas as pd
import geopandas as gpd

def change_wpt(position, capacity, height, turbine):
    args = {
        'lat': 51.8000,  # 51.5000-52.0000
        'lon': 12.2000,  # 11.8000-13.1500
        'date_from': '2011-01-01',
        'date_to': '2011-12-31',
        'capacity': 1000.0,
        'height': 100,
        'turbine': 'Vestas V164 7000',
        'format': 'json',
        'local_time': 'true',
        'raw': 'false',
    }

    args['capacity'] = capacity
    args['height'] = height
    args['lat'] = position[0]
    args['lon'] = position[1]
    args['turbine'] = turbine

    return args

def get_df(args):
    token = 'Please get your own'
    api_base = 'https://www.renewables.ninja/api/'

    s = requests.session()
    # Send token header with each request
    s.headers = {'Authorization': 'Token ' + token}

    url = api_base + 'data/wind'

    r = s.get(url, params=args)

    parsed_response = json.loads(r.text)
    df = pd.read_json(
    json.dumps(parsed_response['data']),orient='index')
    metadata = parsed_response['metadata']
    return df

enercon_production = get_df(change_wpt(
    position,
    capacity=1,
    height=df[["hub_height"]].mean(),
    turbine=enercon)
)

vestas_production = get_df(change_wpt(
    position,
    capacity=1000,
    height=df[["hub_height"]].mean(),
    turbine=vestas)
)
```

#### Gewichtung und Skalierung der Zeitreihen

Um die Charakteristika der beiden o.g. Anlagentypen zu berücksichtigen, erfolgt
eine gewichtete Summierung der Zeitreihen anhand der berechneten Häufigkeit.

#### Zukunftsszenarien

Analog zu dem oben beschriebenen Vorgehen wird eine separate Zeitreihe für
zukünftige WEA berechnet. Hierbei wird eine Enercon E126 6500 mit einer
Nabenhöhe von 159 m angenommen
([PV- und Windflächenrechner](https://zenodo.org/record/6794558)).

Da die Zeitreihe sich nur marginal von der obigen Status-quo-Zeitreihe
unterscheidet, wird letztere sowohl für den Status quo als auch die
Zukunftsszenarien verwendet.

- Einspeisezeitreihe: `wind_feedin_timeseries.csv`

### Freiflächen-Photovoltaik

#### PV-Anlage (2022)

Stündlich aufgelöste Zeitreihe der Photovoltaikeinspeisung über 1 Jahr auf Basis
von [MaStR](../../digipipe/store/raw/bnetza_mastr/dataset.md) und
[renewables.ninja](http://renewables.ninja).
Wie bei der Windeinspeisung wird auf eine Auflsöung auf Gemeindeebene aufgrund
geringer regionaler Abweichungen verzichtet.

Für die Generierung der Zeitreihe über
[renewables.ninja](http://renewables.ninja)
wird eine Position(lat, lon), Nennleistung (capacity), Verluste (system_loss)
Nachführung (tracking), Neigung (tilt) und der Azimutwinkel (azim) benötigt.

Als Position wird analog zur Windenergieanlage der räumlicher Mittelwert
verwendet. Laut MaStR werden lediglich 13 Anlagen nachgeführt (0,01 % der
Kapazität), die Nachführung wird daher vernachlässigt. Die Neigung ist aus MaStR
nicht bekannt, es dominieren jedoch Anlagen auf Freiflächen sowie Flachdächern
im landwirtschaftlichen Kontext. Nach
[Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)
wird diese mit 30° angenommen.
Die Nennleistung Wird auf 1 MW gesetzt/normiert.

#### Zukunftsszenarien

Die Status-quo-Zeitreihe wird sowohl für den Status quo als auch die
Zukunftsszenarien verwendet.

- Einspeisezeitreihe: `pv_feedin_timeseries.csv`

### Solarthermie

- Einspeisezeitreihe: `st_feedin_timeseries.csv` (Kopie von
  PV-Einspeisezeitreihe)

### Laufwasserkraft

Hier wird eine konstante Einspeisung angenommen.

- Einspeisezeitreihe: `ror_feedin_timeseries.csv`

**Dataset: `raw/renewables.ninja_feedin`**

??? metadata "Metadata"
    ```json
    {
        "name": "renewables.ninja_feedin",
        "title": "EE-Einspeisezeitreihen",
        "id": "renewables.ninja_feedin",
        "description": "Einspeisezeitreihen f\u00fcr Erneuerbare Energien, normiert auf 1 MW bzw. 1 p.u. Als Wetterjahr wird 2011 verwendet",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Erneuerbare",
            "Energien",
            "Einspeisezeitreihen",
            "renewables.ninja"
        ],
        "publicationDate": "2016-09-21",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2023-04-14",
            "timeseries": [
                {
                    "start": "2011-01-01T00:00+01",
                    "end": "2011-12-31T23:00+01",
                    "resolution": "1 h",
                    "alignment": "left"
                }
            ]
        },
        "sources": [
            {
                "title": "renewables.ninja_feedin",
                "description": "Einspeisezeitreihen f\u00fcr Erneuerbare Energien, normiert auf 1 MW bzw. 1 p.u. Als Wetterjahr wird 2011 verwendet",
                "path": "hhttps://www.renewables.ninja/",
                "licenses": [
                    {
                        "name": "CC BY-NC 4.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://creativecommons.org/licenses/by-nc/4.0/",
                        "instruction": "you are free to copy, redistribute and adapt them for non-commercial purposes, provided you give appropriate credit. Note that the data is made available as-is and without warranty. We cannot guarantee its accuracy, and accept no responsibility for any liability arising from its use. You are advised to examine the quality of the data for your intended purposes, and to consult the publications linked on this page.",
                        "attribution": "\u00a9 www.renewables.ninja, 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "CC BY-NC 4.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://creativecommons.org/licenses/by-nc/4.0/",
                "instruction": "you are free to copy, redistribute and adapt them for non-commercial purposes, provided you give appropriate credit. Note that the data is made available as-is and without warranty. We cannot guarantee its accuracy, and accept no responsibility for any liability arising from its use. You are advised to examine the quality of the data for your intended purposes, and to consult the publications linked on this page.",
                "attribution": "\u00a9 www.renewables.ninja, 2023"
            }
        ],
        "contributors": [
            {
                "title": "aaronschilling",
                "email": "Aaron.Schilling@rl-institut.de",
                "date": "2023-09-07",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Verwaltungsgebiete Deutschlands

Verwaltungsgebiete Deutschlands (Verwaltungsgebiete 1:250 000).

**Dataset: `raw/bkg_vg250`**

??? metadata "Metadata"
    ```json
    {
        "name": "bkg_vg250",
        "title": "Adminstrative areas of Germany",
        "id": "bkg_vb250",
        "description": "Geopackage with administative areas of Germany - Verwaltungsgebiete 1:250 000",
        "language": [
            "en-GB",
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "adminstrative areas",
            "Verwaltungsgebiete"
        ],
        "publicationDate": "2022-01-01",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": "1:250 000"
        },
        "temporal": {
            "referenceDate": "2022-01-01",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Bundesamt f\u00fcr Kartographie und Geod\u00e4sie - Verwaltungsgebiete 1:250 000 VG250 (Ebenen)",
                "description": "Dieser Datensatz stellt die Verwaltungsgebiete 1:250 000 (VG250) mit Stand 01.01. f\u00fcr das Gebiet der Bundesrepublik Deutschland bereit.",
                "path": "https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-250-000-stand-01-01-vg250-01-01.html",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0",
                        "path": "http://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": " \u00a9 GeoBasis-DE / BKG - 2022"
                    }
                ]
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-23",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```
