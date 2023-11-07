# Bevölkerungsentwicklung

EinwohnerInnen je Gemeinde: Historische Daten und Prognosen

## Historische Daten bis 2022

Statistisches Bundesamt (Raw dataset:
[destatis_gv](../../raw/destatis_gv/dataset.md))

## Prognosen bis 2035

Statistisches Landesamt Sachsen-Anhalt (Raw dataset:
[stala_st_pop_prog](../../raw/stala_st_pop_prog/dataset.md)). Deaktivieren
mittels entfernen der Zieljahre in [config.yml](config.yml) im Abschnitt
`prognosis_fstate_munlevel`.

Kann für andere Regionen auch durch DemandRegio (s.u.) ersetzt werden, die
tatsächliche regionale Auflösung wird dadurch reduziert.

## Prognosen bis 2045

DemandRegio (Raw dataset: [demandregio](../../raw/demandregio/dataset.md))
basierend auf der
[14. koordinierten Bevölkerungsvorausberechnung](https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Bevoelkerungsvorausberechnung/aktualisierung-bevoelkerungsvorausberechnung.html)
der Statistischen Ämter von Bund und Ländern. Diese Daten liegen auf
Landkreisebene vor, daher erfolgt eine gleichmäßige Skalierung der
dazugehörigen Gemeinden auf den jeweiligen Prognosewert.

Deaktivieren mittels entfernen der Zieljahre in [config.yml](config.yml) im
Abschnitt `prognosis_germany_districtlevel`.

## Extrapolation

Über 2045 hinaus wird lineare Extrapolation auf Basis der letzten beiden
Prognosejahre unterstützt. Um diese zu aktivieren, müssen lediglich Zieljahre
in die [config.yml](config.yml) im Abschnitt `extrapolation` eingetragen werden.
