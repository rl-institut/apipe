# Potenzialgebiete PV-Freiflächen (Deutschland gesamt)

Diese Analyse widmet sich der Identifizierung von Photovoltaik-Potenzialen auf
Agrarflächen, basierend auf einer Auswertung rasterbasierter Landnutzungsdaten.
Im Fokus stehen dabei die Bodenqualität und die Eignung landwirtschaftlicher
Dauerkulturen für die Integration von Agri-PV-Systemen.

## Datengrundlage und Methodik

Datengrundlage:

- SQR-Daten (Soil Quality Rating), Datensatz
  [bgr_sqr](../../raw/bgr_sqr/dataset.md).
  - **(I)** SQR Originaldaten
- Potenzialflächen Agri-PV, Datensatz
  [oei_agri_pv](../../raw/oei_agri_pv/dataset.md)
  - **(II)** SQR Gesamtpotenzial (`Agri-PV-Potenziale_Gesamt_100x100_EPSG3035`)
  - **(III)** SQR 50-70
    (`Agri-PV-Potenziale_SQR_50-70_100x100_EPSG3035`)
- Feldblockkataster Brandenburg, Datensatz
  [mluk_bb_field_block_cadastre](../../preprocessed/mluk_bb_field_block_cadastre/dataset.md)
  - **(IV)** Dauerkulturen (`DFBK_FB.tif`)

(II) und (IV) sowie (III) und (IV) sind nicht disjunkt, für die Erstellung von
**(A)** und **(B)** erfolgt eine Differenzbildung, s.u.

Es werden folgende Flächen verwendet, auf welchen jeweils eine andere
technologische Umsetzung als Grundlage angenommen wird:

- **(A)** Auf Acker- und Grünlandflächen mit sehr geringer..geringer Bodengüte
  (SQR 0..50): Klassische, niedrig aufgeständerte FF-PV.
- **(B)** Auf Acker- und Grünlandflächen mit geringer..mittlerer Bodengüte (SQR
  50..70): Agri-PV - bifaziale, vertikal aufgeständerte PV
- **(C)** Dauerkulturen: Agri-PV - hoch aufgeständerte PV-Systeme

Berechnung/Verschneidung

- **(A)** = **(II)**-(**(I)** mit Wert>50)-(**(IV)** mit Wert>0)
- **(B)** = **(III)**-(**(IV)** mit Wert>0)
- **(C)** = **(IV)**

## Ergebnisdaten

- **(A)**: `potentialarea_pv_ground_soil_quality_low.tif`
- **(B)**: `potentialarea_pv_ground_soil_quality_medium.tif`
- **(C)**: `potentialarea_pv_ground_permanent_crops.tif`
