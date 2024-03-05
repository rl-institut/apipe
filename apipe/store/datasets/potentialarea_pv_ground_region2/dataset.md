# Regionale Analyse von Agri-PV-Potenzialen

Diese Analyse widmet sich der Identifizierung von Photovoltaik-Potenzialen auf
Agrarflächen, basierend auf einer Auswertung rasterbasierter Landnutzungsdaten.
Im Fokus stehen dabei die Bodenqualität und die Eignung landwirtschaftlicher
Dauerkulturen für die Integration von Agri-PV-Systemen.

## Methodik und Datengrundlage

Grundlage der Analyse bilden SQR-Daten (Soil Quality Rating) des Bundesinstituts
für Geowissenschaften und Rohstoffe (BGR), ergänzt durch spezifische
Rasterdaten, die Agri-PV-Potenziale illustrieren. Die Datenkategorien umfassen:

- **(A) SQR Originaldaten**: Allgemeine Bodenqualität, basierend auf dem
  Müncheberger Soil Quality Rating (Müller et al., 2008).
- **(B) Gesamte Agri-PV-Potenzialfläche**: Flächen mit Gesamtpotenzial für
  Agri-PV.
- **(C) Agri-PV-Potenzialfläche (SQR 50-70)**: Fokus auf Böden mit geringer bis
  mittlerer Bodengüte für bifaziale, vertikal aufgeständerte PV-Anlagen.
- **(D) Dauerkulturen**: Potenzialflächen für hochaufgeständerte PV-Systeme.

Die Methodik schließt die Verschneidung und Auswertung der Rasterdaten ein, um
spezifische PV-Potenzialbereiche zu identifizieren, die nachfolgend vektorisiert
und im Hinblick auf Mindestflächengröße sowie durchschnittliche Rasterwerte pro
Polygon analysiert werden.

**Quellen**:

- Müller, L., et al. (2008). Das Müncheberger Soil Quality Rating (SQR).
  Berichte der DBG.
- BGR (2023). Nutzungsdifferenzierte Bodenübersichtskarte von Deutschland 1:
  1.000.000. [BGR Bodenübersichtskarte](https://www.bgr.bund.de/DE/Themen/Boden/Ressourcenbewertung/Ertragspotential/Ertragspotential_node.html).

## Ergebnisse

Die Ergebnisse bieten Einblicke in das regionale PV-Potenzial und umfassen:

- **GeoPackages**: Vektorisierte Darstellung der Agri-PV-Potenzialflächen,
  differenziert nach Bodenqualität und Kulturtyp:
    - `potentialarea_pv_ground_soil_quality_low_region.gpkg`
    - `potentialarea_pv_ground_soil_quality_medium_region.gpkg`
    - `potentialarea_pv_ground_permanent_crops_region.gpkg`

### Statistische Auswertung

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet, siehe `potentialarea_pv_ground_area_stats_muns.csv`. Die
Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

### Ausbauziele

Es werden regionalisierte PV-Ausbauziele für die Region berechnet, indem die
Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der regional verfügbaren Potenzialflächen disaggregiert werden. Hierzu
wird der Anteil der Flächensumme der drei o.g. Flächentypen an den bundesweit
verfügbaren Flächen (Datensatz [oei_agri_pv](../../raw/oei_agri_pv/dataset.md))
berechnet. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen- und
Aufdach-PV unterschieden wird, wird ein Verhältnis von 50:50 angenommen, d.h.
es entfallen bundesweit 214 GW auf Freiflächen-PV (inkl. Agri-PV).

Ergebnisse in `potentialarea_pv_ground_regionalized_targets.json`
