# Emissionen
Emissionen für die Jahre 1990 und 2019 für Sachsen-Anhalt (aus
[THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf)) und disaggregiert für die Region ABW.

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

## Disaggregation
Anhand unterschiedlicher Kriterien und Datenquellen wurde näherungsweise von den vorliegenden Emissionen für Sachsen-Anhalt für 1990 und 2019 auf die Region ABW disaggregiert. Je Sektor sind hier die gewählten **energiebestimmenden Größen (EnbG)** angegeben, sowie die Herangehensweise zur jeweiligen Berechnung.

### Sektor Energiewirtschaft (CRF 1.A.1 + 1.B)
Aus der Liste der [Emissionshandelspflichtigen Anlagen](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3) wurden jene Daten zu Anlagen extrahiert, welche sich in Sachsen-Anhalt befinden und als Bezeichnung "Energieumwandlung >= 50 MW FWL" oder "Energieumwandlung 20–50 MW FWL" (Haupttätigkeit nach TEHG) aufweisen.
Die Summe der angegebenen Emissionen (t CO2 Äq) jener Anlagen, welche in der Region ABW liegen, wurde in Relation zu der Summe der Emissionen aus den Anlagen in Gesamt-ST gesetzt. Dieser Anteil wurde auf die im THG-Bericht angegebene Emissionsmenge im Sektor "Energiewirtschaft (1.A.1)" sowie "Prozessemissionen (1.B)" angelegt und so für ABW näherungsweise disaggregiert.

Hinweise:
- Aufgrund mangelnder Daten wurde für das Jahr 1990 auf die ältesten verfügbaren Daten (2005-2007) aus der Anlagenliste zurückgegriffen.
- Energiewirtschaftlich relevante Anlagen unter 20 MW FWL sind in der Anlagenliste nicht erfasst und konnten somit nicht berücksichtigt werden.

Quellen:
- [Emissionshandelspflichtige Anlagen in Deutschland 2020 (Stand 03.05.2021)](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3)
- [Treibhausgasemissionen in Sachsen-Anhalt 2018 (Stand 12.05.2021)](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/THG_Bericht_2018.pdf)

##### CRF 1.A.1
Energiewirtschaft (Umwandlungsbereich): umfasst die öffentliche Elektrizitäts- und Wärmeversorgung sowie Raffinerien.

EnbG: Emissionen aus europäischem Emissionshandel

##### CRF 1.B
Diffuse Emissionen aus Brennstoffen: Diese Kategorie beinhaltet flüchtige Emissionen aus der Gewinnung, Verarbeitung und Verteilung von Brennstoffen. Die wichtigsten Quellen sind die Verteilung von Erdgas, aber auch Emissionen aus Förderung und Abfackelung, die Extraktion und Umwandlung von Braunkohle, Emissionen aus der Raffination von Erdöl sowie Emissionen aus der Lagerung und Verteilung von Mineralölprodukten.

EnbG: Emissionen aus europäischem Emissionshandel

### Sektor Industrie (CRF 1.A.2)
Dieser Sektor umfasst sämtliche energiebedingten Emissionen durch Verarbeitendes Gewerbe.

Zur Disaggregierung wurde der Energieverbrauch der Industriebetriebe in ABW mit dem Gesamtenergieverbrauch aller Industriebetriebe in Sachsen-Anhalt in Relation gesetzt.
Dabei wurde eine Differenzierung hinsichtlich der Energieträgerzusammensetzung von ABW im Vergleich zu ST durchgeführt und anhand etablierter Emissionsfaktoren berechnet.

EnbG: Energieverbrauch nach Energieträgern

Quellen:
- [Energieverbrauch der Industriebetriebe in Sachsen-Anhalt nach ausgewählten Energieträgern und Kreisen](https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Energie/Tabellen/Energieverwendung/Energieverbrauch_nach_Kreisen_ab_dem_Jahr_2010.xlsx)
- [Emissionsfaktor für Stromerzeugung (UBA)](https://www.umweltbundesamt.de/sites/default/files/medien/479/bilder/dateien/entwicklung_der_spezifischen_emissionen_des_deutschen_strommix_1990-2020_und_erste_schaetzungen_2021.pdf)
- [BISKO Bilanzierungs-Systematik Kommunal (Aktualisierung 11/2019)](https://www.ifeu.de/fileadmin/uploads/BISKO_Methodenpapier_kurz_ifeu_Nov19.pdf)

### Sektor Prozessemissionen (CRF 2)
Dieser Sektor umfasst sämtliche Emissionen, welche durch diverse Industriprozesse anfallen. Konkreter sind das Emissionen aus: Herstellung mineralischer Produkte, chemischer Industrie, Herstellung von Metallen, übrigen Prozessen und Produktverwendungen. (CRF 2.A-H)
Zur Disaggregierung wurde erneut die [Liste der Emissionshandelspflichtigen Anlagen](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3) herangezogen. In diesem Fall wurde allerdings der Anteil aller Anlagen, welche nicht der Energiewirtschaft zugerechnet werden, zur Bestimmung des Anteils von AWB an ST gewählt.

EnbG: Emissionen aus europäischem Emissionshandel

### Sektor Verkehr (CRF 1.A.3)
Dieser Sekotr umfasst Emissionen aus dem Straßenverkehr, dem zivilen Luftverkehr, aus dem Schiffsverkehr, verbrennungsbedingte Emissionen aus dem Schienenverkehr sowie Emissionen des übrigen Verkehrs und weitere Quellen zur Bereitstellung der im Verkehr verbrauchten Energie.
Die Verbrennung von Mineralölprodukten im Straßenverkehr spielt die größte Rolle und macht weit über 90 % der sektoralen Emissionen aus.
Daher wird zur Disaggreagation der motorisierte Straßenverkehr über zugelassene Kraftfahrzeuge mit durchschnittlichen Fahrleistungen und spezifischer Emissionen pro Kilometer und Fahrzeugklasse herangezogen.

Zunächst wird aus [Verkehr in Kilometern (VK) ZeitreiheJahre 2014 - 2022](https://www.kba.de/DE/Statistik/Kraftverkehr/VerkehrKilometer/vk_inlaenderfahrleistung/vk_inlaenderfahrleistung_node.html;jsessionid=DD419FD0604C0BCC72A9E4533BB0319F.live21324) und [Umweltfreundlich mobil! Ein ökologischer Verkehrsartenvergleich für den Personen- und Güterverkehr in Deutschland)](https://www.umweltbundesamt.de/sites/default/files/medien/5750/publikationen/2021_fb_umweltfreundlich_mobil_bf.pdf) ein durchschnittlicher CO2-Emissionswert pro Jahr und Fahrzeugklasse ermittelt. Dieser wird dann mit den zugelassenen Fahrzeugen der entsprechenden Fahrzeugklassen aus [Kraftfahrzeugbestand nach Kraftfahrzeugarten - Stichtag 01.01. - regionale Tiefe: Kreise und krfr. Städte (bis 01.01.2019)](https://www-genesis.destatis.de/genesis//online?operation=table&code=46251-0001&bypass=true&levelindex=0&levelid=1691405772899#abreadcrumb) einmal für ganz Sachsen-Anhalt und einmal ABW multipliziert. Daraus kann dann ein Verhältnis gewonnen werden, dass den prozentualen Anteil der Verkehrsemissionen von ABW kommt.
Dieser prozentuale Anteil wird auf die Verkehrsemissionen aus dem [THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf) angewednet.

Hinweise:
- Die Datenlage für die zugelassenen Fahrzeuge, gefahrenen Kilometer und Emissionen pro km sind nicht spezifisch für 1990 sondern nur für einzelne Jahre der frühen 1990er verfügbar. Daher ist der Emissionswert für 1990 it einer höheren Unsicherheit behaftet.

EnbG:
* Zugelassene Kraftfahrzeuge
* gewichtet mit durchschn. Fahrleistung und spez. CO2 Emission pro km und Fahrzeugklasse

Quellen:
- [Kraftfahrzeugbestand nach Kraftfahrzeugarten - Stichtag 01.01. - regionale Tiefe: Kreise und krfr. Städte (bis 01.01.2019)](https://www-genesis.destatis.de/genesis//online?operation=table&code=46251-0001&bypass=true&levelindex=0&levelid=1691405772899#abreadcrumb)
- [Umweltfreundlich mobil! Ein ökologischer Verkehrsartenvergleich für den Personen- und Güterverkehr in Deutschland)](https://www.umweltbundesamt.de/sites/default/files/medien/5750/publikationen/2021_fb_umweltfreundlich_mobil_bf.pdf)
- [Verkehr in Kilometern (VK) ZeitreiheJahre 2014 - 2022](https://www.kba.de/DE/Statistik/Kraftverkehr/VerkehrKilometer/vk_inlaenderfahrleistung/vk_inlaenderfahrleistung_node.html;jsessionid=DD419FD0604C0BCC72A9E4533BB0319F.live21324)

### Sektor Sonstige Energie (insbes. Gebäude) (CRF 1.A.4 + 1.A.5)

Dieser Sektor umfasst den durch Energieumwaldnung nicht bereits abgedeckten Energiebedarf. Das sind vor allem Feuerungsanlagen von kleinen Einzelraumfeuerungen (z. B. Kaminöfen) bis hin zu immissionsschutzrechtlich genehmigungsbedürftigen Anlagen mit einer Nennwärmeleistung von mehreren Megawatt.
Zur Disaggreagtion wurde daher der Wärmebedarf von ABW im Verhältnis zum Wärmebedarf von gesamt Sachsen Anhalt gewählt. Der Wärmevedarf umfasst Raumwärme, Warmwasser und Kochen und wird aus Daten aus dem  Pipeline-Datensatz demand_heat_region generiert, mehr Details siehe dort.

Ergebnis: 17,46 % des Bedarfs in Sachsen-Anhalt entfällt auf ABW.
Dieser Wert wird mit den Emissionen auf den [THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf) angewednet.

Code
```
# Sektor HH
heat_hh_dist_states = gpd.read_file("demand_heat_zonal_stats-res-bkg_vg250_federal_states.gpkg")
heat_hh_demand_st = float(heat_hh_dist_states.loc[heat_hh_dist_states.nuts == "DEE"].heat_demand)
heat_hh_demand_abw = gpd.read_file("demand_heat_zonal_stats-res-bkg_vg250_muns_region.gpkg").heat_demand.sum()

# Sektor GHD
heat_cts_dist_states = gpd.read_file("demand_heat_zonal_stats-ser-bkg_vg250_federal_states.gpkg")
heat_cts_demand_st = float(heat_cts_dist_states.loc[heat_cts_dist_states.nuts == "DEE"].heat_demand)
heat_cts_demand_abw = gpd.read_file("demand_heat_zonal_stats-ser-bkg_vg250_muns_region.gpkg").heat_demand.sum()

# Anteil ABW an ST
heat_share = (heat_hh_demand_abw + heat_cts_demand_abw) / (heat_hh_demand_st + heat_cts_demand_st)
```

EnbG: Wärmebedarf aus Energiesystem


### Sektor Landwirtschaft (CRF 3)
Der Sektor umfasst Emissionen aus der Viehwirtschaft und der Bewirtschaftung von Böden.
Daher werden zunächst die Emissionsunterkategorien 3.A-J aus dem [THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf) Viehwirtschaft oder der Bewirtschaftung von Böden zugeordnet. Danach werden diese getrennt nach den Viehbeständen bzw. der landwirtschaftlich genutzen Fläche disaggreiert.


##### CRF 3.A - Landwirtschaft – Fermentation
Emission durch Fermentation (CRF 3.A) entseht durch Verdauungsprozesse in der Viehwirtschaft. Deswegen kann der Anteil ABWs an diesen Emissionen durch die Viehbestände aus [Viehbestand der landwirtschaftlichen Betriebe in Großvieheinheiten (GV) nach Jahren und Kreisen)](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/land-und-forstwirtschaft-fischerei/tabellen-viehwirtschaft-und-tierische-erzeugnisse#c234218) abgeschätzt werden.

Hinweise:
- Die Viehbestände für 1990 sind nicht bekannt, es wird stattdessen auf die Viehbestände von 1996 zurückggegriffen.

EnbG: Viehbestände


Quellen:
- [Viehbestand der landwirtschaftlichen Betriebe in Großvieheinheiten (GV) nach Jahren und Kreisen)](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/land-und-forstwirtschaft-fischerei/tabellen-viehwirtschaft-und-tierische-erzeugnisse#c234218)

##### CRF 3.B-J:
Die Unterkategorien 3.C-J ist eine Proportionalität der Emissionen und der landwirtschafltich genutzen Fläche zu erwarten. Unterkategorie 2.B "Wirtschaftsdüngerausbringung (ohne Gärreste)" ist allerdings ein Grenzfall, da er aus Abfällen der Tierhaltung produziert wird und schon dabei Treibhausgase entstehen, diese aber nicht vor Ort eingesetzt werden müssen, sondern auf beliebigen landwirtschafltichen Flächen eingesetzt werden kann. Daher wird hier auch diese Unterkategorie der Landnutzung zugeordnet.
Die Anteile der landwirtschaftlich genutzen Fläche von ABW in Sachsen-Anhalt sind der Tabelle Flaeche_nach_Kultuarten_nach_Jahren_und_Kreisen](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/land-und-forstwirtschaft-fischerei/tabellen-bodennutzung-und-anbau) zu entnehmen.
Dieses Verhältnis kann auf die Emissionen aus dem [THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf) angewendet werdsen.

Hinweis:
- die Flächenntuzungsadaten gehen nicht bis  1990 zurück, ändern sich über die Jahre aber nur marginal, sodass hier auch nicht von großen Abweichungen auszugehen ist.

EnbG: landwirtschaftlich genutzte Fläche


Quellen:
- [Flaeche_nach_Kultuarten_nach_Jahren_und_Kreisen](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/land-und-forstwirtschaft-fischerei/tabellen-bodennutzung-und-anbau)

### Sektor Abfall und Abwasser (CRF 5)
Dieser Sektor besteht vor allem aus Emissionen aus Abfalldeponien, welche der Zersetzung organischer Materialien in Deponien entstehen.
Es wird angenommen, dass der Abfall aus Produktionsprozessen gegenüber den Abfällen aus Konsum vernachlässigbar sind, weswegen eine Disaggregation auf Grundlage der Bevölkerung von ABW vorgenommen wird.

EnbG: Bevölkerung ABW

Quellen:
- [Bevölkerung nach Geschlecht in den Gemeinden](https://genesis.sachsen-anhalt.de/genesis//online?operation=table&code=12411-0001&bypass=true&levelindex=0&levelid=1691507280245#abreadcrumb)
