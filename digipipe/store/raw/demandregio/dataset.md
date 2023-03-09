# DemandRegio

Regionalisierte Bevölkerungsprognose sowie Strom-, Wärme und Gasbedarf auf
Landkreisebene.

Die Daten wurden abgerufen mit einer
[modifizierten Version des DemandRegio disaggregators](https://github.com/nesnoj/disaggregator), 
in der softwareseitige, jedoch keine methodischen Änderungen vorgenommen wurden.

Der disaggregator basiert auf Daten bis 2017, anschließende Jahre werden
fortgeschrieben.

Weitere Informationen zum Projekt DemandRegio
- [Abschlussbericht](https://www.ffe.de/wp-content/uploads/2020/10/DemandRegio_Abschlussbericht.pdf)
- [Abschlussworkshop](https://www.tu.berlin/er/forschung/projekte/demandregio-2)

Die erzeugten Rohdaten wie unten beschrieben wurden mittels 
[API](http://opendata.ffe.de:4000/) abgerufen. Diese können alternativ direkt
vom [OpenData-Portal der FfE](https://opendata.ffe.de/project/demandregio/)
bezogen werden.

**Installation (digipipe venv aktiviert):**

```commandline
pip install openpyxl==3.1.0
pip install disaggregator@git+https://github.com/nesnoj/disaggregator.git#egg=disaggregator
```

Annahmen und Parameter
- Wetterjahr: Einheitlich 2011

## Details zum Datenabruf

### Bevölkerung

Bevölkerung (Summe) und Bevölkerung je Haushaltsgröße (1, 2, 3, 4, 5, >5) je
NUTS3.

Jahre
- Bevölkerung bis 2017 historische Werte
- Bevölkerung ab 2018 prognostizierte Werte basierend auf der 14. koordinierten
  Bevölkerungsvorausberechnung der Statistischen Ämter von Bund und Ländern.
- Haushalte nur 2011

```python
import pandas as pd
from disaggregator import data

# Population
dr_hh_population = pd.DataFrame()
for year in [2010, 2015, 2017, 2020, 2021, 2022, 2025, 2030, 2035, 2040, 2045]:
    dr_hh_population[year] = round(data.population(year=year)).astype(int)

dr_hh_population.to_csv("dr_hh_population.csv")

# Households
data.households_per_size().to_csv("dr_hh_households_2011.csv")
```

### Haushalte: Strom

Bedarfe und SLP-Zeitreihen je NUTS3 mit Bottom-Up-Methode nach Haushaltsgröße.

Jahre
- 2017: Letzte verfügbare Daten
- 2022: Status quo, Fortschreibung mit Berücksichtigung Demografie und
  Wanderung
- 2035: Fortschreibungsjahr mit Berücksichtigung Demografie und Wanderung
- 2045: Fortschreibungsjahr

```python
from disaggregator import spatial, temporal

for year in [2017, 2022, 2035, 2045]:
  # Consumption
  spatial.disagg_households_power(
      by="households",
      weight_by_income=True,
      year=year,
      scale_by_pop=True,
  ).to_csv(f"dr_hh_power_consumption_{year}.csv")
  
  # Timeseries
  temporal.disagg_temporal_power_housholds_slp(
      use_nuts3code=True,
      by="households",
      weight_by_income=True,
      year=year,
      scale_by_pop=True,
  ).to_csv(f"dr_hh_power_timeseries_{year}.csv")
```

### Haushalte: Wärme/Gas

Bedarfe und SLP-Zeitreihen je NUTS3

TBD

### GHD und Industrie: Strom

Bedarfe und Zeitreihen je NUTS3
- Bedarfe: je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: für alle WZ aggregiert, Einzelprofile basieren je nach WZ
  auf gemessenen oder SLP inkl. Wanderung

Jahre
- 2017: Letzte verfügbare Daten
- 2022: Status quo, Fortschreibung mit Berücksichtigung Beschäftigte und
  Effizienzgewinne
- 2035: Max. Fortschreibungsjahr mit Berücksichtigung Beschäftigte und
  Effizienzgewinne

```python
from disaggregator import spatial, temporal

# CTS
for year in [2017, 2022, 2035]:
  # Consumption
  spatial.disagg_CTS_industry(
      sector='CTS',
      source='power',
      use_nuts3code=True,
      year=year,
  ).to_csv(f"dr_cts_power_consumption_{year}.csv")
  # Timeseries
  temporal.disagg_temporal_power_CTS(
      detailed=False,
      use_nuts3code=True,
      year=year,
  ).to_csv(f"dr_cts_power_timeseries_{year}.csv")

# Industry
for year in [2017, 2022, 2035]:
  # Consumption
  spatial.disagg_CTS_industry(
      sector='industry',
      source='power',
      use_nuts3code=True,
      year=year,
  ).to_csv(f"dr_ind_power_consumption_{year}.csv")
  # Timeseries
  temporal.disagg_temporal_industry(
      source="power",
      detailed=False,
      use_nuts3code=True,
      no_self_gen=False,
      year=year,
  ).to_csv(f"dr_ind_power_timeseries_{year}.csv")
```

### GHD und Industrie: Wärme/Gas

Bedarfe und SLP-Zeitreihen je NUTS3

TBD
