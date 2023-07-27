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
Anhand unterschiedlicher Kriterien und Datenquellen wurde näherungsweise von den vorliegenden Emissionen für Sachsen-Anhalt für 1990 und 2019 auf die Region ABW disaggregiert. Je Sektor sind hier die **energiebestimmenden Größen (EnbG)** angegeben, sowie die Herangehensweise zur jeweiligen Berechnung.

### Sektor Energiewirtschaft (CRF 1.A.1 + 1.B)
Aus der Liste der [Emissionshandelspflichtigen Anlagen](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3) wurden jene Daten zu Anlagen extrahiert, welche sich in Sachsen-Anhalt befinden und als Bezeichnung "Energieumwandlung >= 50 MW FWL" oder "Energieumwandlung 20–50 MW FWL" (Haupttätigkeit nach TEHG) aufweisen.
Die Summe der angegebenen Emissionen (t CO2 Äq) jener Anlagen, welche in der Region ABW liegen, wurde in Relation zu der Summe der Emissionen aus den Anlagen in Gesamt-ST gesetzt. Dieser Anteil wurde auf die im THG-Bericht angegebene Emissionsmenge im Sektor "Energiewirtschaft (1.A.1)" sowie "Prozessemissionen (1.B)" angelegt und so für ABW näherungsweise disaggregiert.

Hinweise:
- Aufgrund mangelnder Daten wurde für das Jahr 1990 auf die ältesten verfügbaren Daten (2005-2007) aus der Anlagenliste zurückgegriffen.
- Energiewirtschaftlich relevante Anlagen unter 20 MW FWL sind in der Anlagenliste nicht erfasst und konnten somit nicht berücksichtigt werden.


###### CRF 1.A.1
Energiewirtschaft (Umwandlungsbereich): umfasst die öffentliche Elektrizitäts- und Wärmeversorgung sowie Raffinerien.

EnbG: Emissionen aus europäischem Emissionshandel


###### CRF 1.B
Diffuse Emissionen aus Brennstoffen: Diese Kategorie beinhaltet flüchtige Emissionen aus der Gewinnung, Verarbeitung und Verteilung von Brennstoffen. Die wichtigsten Quellen sind die Verteilung von Erdgas, aber auch Emissionen aus Förderung und Abfackelung, die Extraktion und Umwandlung von Braunkohle, Emissionen aus der Raffination von Erdöl sowie Emissionen aus der Lagerung und Verteilung von Mineralölprodukten.

EnbG: Emissionen aus europäischem Emissionshandel

Quellen:
- [Emissionshandelspflichtige Anlagen in Deutschland 2020 (Stand 03.05.2021)](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2013-2020/2020.pdf?__blob=publicationFile&v=3)
- [Treibhausgasemissionen in Sachsen-Anhalt 2018 (Stand 12.05.2021)](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/THG_Bericht_2018.pdf)


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
EnbG: Emissionen aus europäischem Emissionshandel

### Sektor Verkehr (CRF 1.A.3)

EnbG:
* Zugelassene Kraftfahrzeuge
* gewichtet mit durchschn. Fahrleistung und spez. CO2 Emission pro km und Fahrzeugklasse

### Sektor Sonstige Energie (insbes. Gebäude) (CRF 1.A.4 + 1.A.5)

EnbG: Wärmebedarf aus Energiesystem

### Sektor Landwirtschaft (CRF 3)

###### CRF 3.A - Landwirtschaft – Fermentation

EnbG: Viehbestände

###### CRF 3.B-J:

EnbG: landwirtschaftlich genutzte Fläche

### Sektor Abfall und Abwasser (CRF 5)

EnbG: Bevölkerung ABW
