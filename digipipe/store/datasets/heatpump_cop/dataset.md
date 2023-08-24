# Wärmepumpen COP

Zeitreihe für die Leistungszahl / Coefficient of performance (COP) für
Wärmepumpen. Berücksichtigt werden Luftwärmepumpen (ASHP) und Erdwärmepumpen
(GSHP). Der COP wird mit Hilfe von Zeitreihen der Umgebungstemperatur (ASHP)
bzw. der Bodentemperatur (GSHP) für jeden Zeitschritt berechnet.

Details zur Berechnungsmethodik können der Dokumentation von
[oemof.thermal](https://oemof-thermal.readthedocs.io/en/latest/compression_heat_pumps_and_chillers.html)
entnommen werden.

Annahmen
- Vorlauftemperatur: 40 °C
- Gütegrad / Quality grade: 0.4 (nach
  [VDE](https://www.energiedialog2050.info/wp-content/uploads/simple-file-list/VDE_ST_ETG_Warmemarkt_RZ-web.pdf))
- Vereisungsverluste bei ASHP: 20 % bei <2 °C

Daraus ergibt sich eine mittlere Jahresarbeitszahl (JAZ) von 3,3 für ASHP und
4,3 für GSHP, die mit typischen Werten für 2019
([AEW](https://static.agora-energiewende.de/fileadmin/Projekte/2022/2022-04_DE_Scaling_up_heat_pumps/A-EW_273_Waermepumpen_WEB.pdf))
übereinstimmen. Für das Zukunftsszenario wird ferner ein Effizienzgewinn durch
technische Weiterentwicklung von 25 % angenommen
[ewi](https://www.ewi.uni-koeln.de/cms/wp-content/uploads/2015/12/2014_06_24_ENDBER_P7570_Energiereferenzprognose-GESAMT-FIN-IA.pdf).

Beide separat erstelle Zeitreihen werden anhand der heutigen Marktdurchdringung
gewichtet und in eine mittlere Zeitreihe für Wärmepumpen überführt. Im Jahr
2022 wurden 87 % ASHP und 13 % GSHP abgesetzt nach
[BWP](https://www.waermepumpe.de/fileadmin/user_upload/waermepumpe/05_Presse/01_Pressemitteilungen/BWP_Branchenstudie_2023_DRUCK.pdf),
über die letzten 10 Jahre beträgt das Verhältnis ca. 80:20.
Für 2045 wird daher ein Anteil von 80 % ASHP und 20 % GSHP angenommen.

Verwendet Datensätze
- [dwd_temperature](../../preprocessed/dwd_temperature/dataset.md)
