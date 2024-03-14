# Solaratlas Brandenburg - Photovoltaikanlagen auf Dachflächen (WFBB)

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Brandenburg der Wirtschaftsförderung Berlin-Brandenburg.

Quelle:
[Solaratlas Brandenburg](https://energieportal-brandenburg.de/cms/inhalte/tools/solaratlas-brandenburg/datenservice)

Hinweise:

- Datensatzbeschreibung:
  `solaratlas_datensatzbeschreibung_eignung_dachflaechen_pv.xlsx`
- Das Attribut `dach_potenzial_ertrag_kwh_a` gibt den potenziellen Ertrag in
  **MWh** an, nicht wie in Datensatzbeschreibung angegeben in kWh
- Referenzwert für die Eignung (100 %) sind 1065 Volllaststunden (berechnet aus
  `dach_potenzial_ertrag_kwh_a` / `dach_potenzial_leistung_kw`  * 1000 /
  `dach_eignung_prozent` * 100)

Beschreibung:

> Mit dieser Analyse sind die verfügbaren Flächen für Photovoltaikanlagen das
> prinzipiell realisierungsfähige Potenzial im Land Brandenburg ermittelt
> worden.
> Dabei sind alle theoretischen Potenziale aufgezeigt, von großen Freiflächen
> bis
> hin auf die Ebene von einzelnen Gebäuden.
>
> Die Solarpotenziale auf Dachflächen wurden mit Hilfe dieser Datensätze
> berechnet:
>
> - Amtliches Liegenschaftskatasterinformationssystem (ALKIS): Enthält sämtliche
    Informationen zu den im Land befindlichen Liegenschaften (Gebäude und
    Flurstücke).
> - LOD2: Informationen zur Beschaffenheit der Dachflächen aller Gebäude im Land
    Brandenburg. Hierbei handelt es sich um ein vereinfachtes 3D-Gebäudemodell,
    bei dem jedem Gebäude eine passende standardisierte Dachform zugeordnet ist.
    Etwaige Dachaufbauten wie Kamine, Antennen oder Dachfenster sind in dem
    Datensatz nicht enthalten.
> - bDOM: Bildbasiertes digitales Oberflächenmodell des Landes Brandenburg, das
    für die Verschattungsanalyse genutzt wurde.
> - Photovoltaic Geographical Information System (PVGIS): System der EU, mit dem
    zur Kalibrierung der Einstrahlungswerte einer Dachfläche für jede Neigungs-
    und Ausrichtungskombination ein Verlustfaktor errechnet wurde.

> Die ermittelten potenziell geeigneten Dachflächen wurden anschließend dem
> ALKIS-Objektkatalog der Gebäudenutzung den Hauptkategorien Öffentliche Zwecke,
> Wohnen, Wirtschaft/Gewerbe und Sonstige zugeordnet.
>
> Kleinstgebäude, die für die Errichtung netzgekoppelter Photovoltaikanlagen in
> der Regel nicht in Frage kommen, wurden ausgeschlossen (Dachfläche bei
> Schrägdächern < 3 m², bei Flachdächern < 6 m²). Zur Berechnung der möglichen
> Leistung und Energiemenge wurde ein Referenzmodul gesetzt (300 Wp-Solarmodul mit
> den Abmessungen 1,65 m x 1,0 m, Modulwirkungsgrad von 18 %).
>
> Die geeigneten Dachflächen wurden in Eignungsklassen eingruppiert. Hat ein
> Gebäude mehrere geeignete Dachflächen, so richtet sich die Eignungsklasse des
> gesamten Gebäudes nach der Eignung der größten Fläche. Die Eignungsklassen
> unterscheiden sich nach der nutzbaren Jahreseinstrahlung, die sich aus Neigung,
> Ausrichtung und Verschattung ergibt:
>
> - gut geeignet: 100 % - 80 %
> - geeignet: 80 % - 60 %
> - bedingt geeignet: 60 % - 40 %
> - nicht geeignet: < 40%

> Bautechnische Gegebenheiten wie der Zustand und die Statik der Gebäude- oder
> Denkmalschutzauflagen sowie Einschränkungen durch kommunale Satzungen sind bei
> der Analyse nicht betrachtet worden. Die Ergebnisse der Dachflächen sind somit
> als erste Einschätzung geeignet und ersetzen für die konkrete Anlage nicht die
> Beurteilung durch Fachunternehmen bzw. eine detaillierte Einzelfallprüfung.
