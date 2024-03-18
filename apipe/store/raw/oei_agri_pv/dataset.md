# Agri-PV-Potenzialflächen

Potenzialflächen für Agri-Photovoltaik des Öko-Institut.

## GIS-Info

Die drei TIFF-Dateien enthalten Informationen darüber, wie viel Potenzialfläche
pro Pixel mit einer Größe von einem Hektar vorhanden ist (1 = 1 ha/ha).
Dies wurde erreicht, indem die Positivflächen aus dem Landnutzungsdatensatz mit
den Ausschlussflächen in QGIS überlagert wurden.
Die Flächen wurden ursprünglich in einem Raster mit einer Größe von 20x20 Metern
erfasst und dann zu 1 Hektar großen Flächen zusammengefasst.

## Positivgebiete

- **Raster_gesamt**:
  - Preidl_Landuse (Landnutzungsklassifikation)
    - Klassen: 1 - 19
- **Raster_besonders_geeignete_Kulturen**
  - Preidl_Landuse
    - Klassen: 9 (legumes) ,11 (leeks), 14 (berries), 15 (stonefruits)
- **Raster_geringe_mittlere_Bodengüte**
  - SQR_Germany:
    - Werte: 50-70

Update: Aus diesem Datensatz werden aufgrund der Datengüte keine Dauerkulturen
("Raster_besonders_geeignete_Kulturen") verwendet.

## Ausschlusskritieren

| Kategorie                | Bezeichnung                         | Quelle / Begründung                                                                                                                         |
|--------------------------|-------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| Biosphärengebiete        | Kernzonen                           | BNatSchG § 25: Nur Maßnahmen, die dem Schutz, der Pflege oder der Entwicklung dienen, sind erlaubt.                                         |
| Biotope                  | Geschützte Biotope                  | BNatSchG § 30: Eingriffe, die die Funktion beeinträchtigen, sind verboten.                                                                  |
| Flughäfen und Flugplätze | Flächen für Flugverkehr             | Luftverkehrsgesetz (LuftVG) § 8: Anlagen dürfen den Flugverkehr nicht gefährden.                                                            |
| Geländeneigung           | Hangneigung größer 20 %             | Keine spezifischen gesetzlichen Beschränkungen, jedoch erhöhte technische Anforderungen und Risiken.                                        |
| Gewässer                 | Fließgewässer und stehende Gewässer | Wasserhaushaltsgesetz (WHG) § 35: Anlagen dürfen Gewässer und deren Nutzung nicht beeinträchtigen.                                          |
| Landschaftsschutzgebiete |                                     | BNatSchG § 26: Anlagen dürfen den besonderen Schutzzweck nicht beeinträchtigen. Aber prinzipiell möglich.                                   |
| Nationale Naturmonumente |                                     | BNatSchG § 22: Anlagen könnten verboten sein, wenn sie den Schutzzweck des Naturmonuments beeinträchtigen.                                  |
| Nationalpark             | Nationalpark                        | BNatSchG § 24: Eingriffe, die die Ziele gefährden, sind verboten.                                                                           |
| Natura 2000-Gebiete      | FFH & Vogelschutzgebiete            | Errichtung nach Baurecht im Einzelfall möglich, allerdings nach EEG § 37 Abs. 3. Aus der Förderkulisse ausgeschlossen.                      |
| Naturdenkmale            | Flächenhafte Naturdenkmale          | BNatSchG § 28: Eingriffe, die die Schutzzwecke beeinträchtigen, sind verboten.                                                              |
| Naturschutzgebiete       | Naturschutzgebiete                  | Bundesnaturschutzgesetz (BNatSchG) § 23: Anlagen dürfen Erhaltungsziele oder den Schutzzweck nicht gefährden.                               |
| Schienenstrecken         | Bahnstrecken                        | Eisenbahn-Bau- und Betriebsordnung (EBO) § 4: Anlagen dürfen den Bahnbetrieb nicht gefährden oder behindern.                                |
| Schienenstrecken         | Bahnverkehrsanlagen                 | Siehe Bahnstrecken.                                                                                                                         |
| Siedlungsflächen         | Wohngebäude                         | Laut Bundesimmissionsschutzgesetz (BImSchG) § 3 sind Schutzabstände zu Wohngebieten notwendig, um die Immissionen zu minimieren.            |
| Siedlungsflächen         | Industrie und Gewerbegebiete        | Keine spezifischen gesetzlichen Beschränkungen, jedoch nicht vorrangige Gebietskulisse für APV. Außerdem nicht in Flächenkulisse enthalten. |
| Straßen                  | Bundesautobahnen                    | Beibehaltung des 15m-Korridors nach EEG 2021 § 37 Abs. 1                                                                                    |
| Straßen                  | Weitere Straßen                     | Siehe Bundesautobahnen.                                                                                                                     |
| Straßen                  | Wege                                | Nicht in Gebietskulisse enthalten                                                                                                           |
| Überschwemmungsgebiete   | Überflutungsflächen HQ100           | WHG § 76: Anlagen dürfen das Hochwasserrisiko nicht erhöhen.                                                                                |
| Wald- und Forstflächen   | Wald                                | Bundeswaldgesetz (BWaldG) § 9: Umwandlung von Wald in andere Nutzungsarten ist für Photovoltaikanlagen verboten.                            |
| Wasserschutzgebietszonen | Zone I und II                       | WHG § 51: Anlagen dürfen die Wasserqualität nicht gefährden.                                                                                |

Datenherkunft:

| Kategorie                                         | Datenquelle                                                                                                                                                   |
|---------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Biosphärengebiete                                 | BfN                                                                                                                                                           |
| Biotope                                           | BfN                                                                                                                                                           |
| Biotopverbund                                     | BfN                                                                                                                                                           |
| Flächenhafte Naturdenkmale                        | BfN                                                                                                                                                           |
| Geländeneigung                                    | EU-DGM                                                                                                                                                        |
| Gewässer                                          | DLM250                                                                                                                                                        |
| Landschaftsschutzgebiete                          | BfN                                                                                                                                                           |
| Landwirtschaftliche Schläge mit Nutzungskategorie | Preidl, Sebastian; Lange, Maximilian; Doktor, Daniel (2020): Land cover classification map of Germany's agricultural area based on Sentinel-2A data from 2016 |
| Nationale Naturmonumente                          | BfN                                                                                                                                                           |
| Nationalpark                                      | BfN                                                                                                                                                           |
| Natura 2000-Gebiete                               | BfN                                                                                                                                                           |
| Naturparke                                        | BfN                                                                                                                                                           |
| Naturschutzgebiete                                | BfN                                                                                                                                                           |
| Sentinel-2A Landnutzungsklassifikation            | Preidl, Sebastian; Lange, Maximilian; Doktor, Daniel (2020): Land cover classification map of Germany's agricultural area based on Sentinel-2A data from 2016 |
| Siedlungsflächen                                  | DLM250                                                                                                                                                        |
| Soil Quality Rating                               | BGR                                                                                                                                                           |

## Quelle

TODO: INSERT ÖI ZENODO SOURCE
