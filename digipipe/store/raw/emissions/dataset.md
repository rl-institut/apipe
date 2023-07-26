# Emissionen
Emissionen für die Jahre 1990 und 2019 für Sachsen-Anhalt (aus
[THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf))
und disaggregiert für die Region ABW.

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
Anhand unterschiedlicher Kriterien und Datenquellen wurde näherungsweise von den vorliegenden Emissionen für Sachsen-Anhalt für 1990 und 2019 auf die Region ABW disaggregiert.

### Sektor Energiewirtschaft (CRF 1.A.1 + 1.B)
Aus der Liste der Emissionshandelspflichtigen Anlagen wurden jene Daten zu Anlagen extrahiert, welche sich in Sachsen-Anhalt befinden und als Angabe "Energieumwandlung >= 50 MW FWL" (Haupttätigkeit nach TEHG) haben.
Die Summe der angegebenen Emissionen (t CO2 Äq) jener Anlagen, welche in der Region AWB liegen, wurde in Relation zu der Summe der Emissionen aus den Anlagen in Gesamt-ST gesetzt. Dieser Anteil wurde auf die im THG-Bericht angegebene Emissionsmenge im Sektor "Energiewirtschaft" angelegt und so für ABW näherungsweise disaggregiert.

Hinweis: Für das Jahr 1990 wurde aufgrund mangelnder Daten auf die ältesten verfügbaren Daten (2005-2007) aus der Anlagenliste zurückgegriffen.

Quelle: [Emissionshandelspflichtigen Anlagen in Deutschland (bis 2020)](https://www.dehst.de/SharedDocs/downloads/DE/anlagenlisten/2020.pdf?__blob=publicationFile&v=3)

###### CRF 1.A.1
EnbG: Emissionen aus europäischem Emissionshandel

###### CRF 1.B
EnbG: Emissionen aus europäischem Emissionshandel

### Sektor Industrie (CRF 1.A.2 + 2)
Hierfür wurde der Energieverbrauch der Industriebetriebe in AWB mit dem Gesamtenergieverbrauch aller Industriebetriebe in Sachsen-Anhalt in Relation gesetzt.

Quellen:
- [Emissionsfaktor für Stromerzeugung (UBA)](https://www.umweltbundesamt.de/sites/default/files/medien/479/bilder/dateien/entwicklung_der_spezifischen_emissionen_des_deutschen_strommix_1990-2020_und_erste_schaetzungen_2021.pdf)
- [Energieverbrauch der Industriebetriebe in Sachsen-Anhalt nach ausgewählten Energieträgern und Kreisen](https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Energie/Tabellen/Energieverwendung/Energieverbrauch_nach_Kreisen_ab_dem_Jahr_2010.xlsx)

###### CRF 1.A.2
EnbG: Energienutzung nach Energieträgern

###### CRF 2 Prozessemissionen
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
