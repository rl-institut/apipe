# DemandRegio

Regionalisierte Bevölkerungsprognose, Haushalte sowie Strom- und Gasbedarfe
inkl. Zeitreihen auf Landkreisebene.

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

Verwendetes Wetterjahr für Gasbedarfszeitreihen: 2011

**Installation (in separater venv):**

```commandline
pip install disaggregator@git+https://github.com/nesnoj/disaggregator.git#egg=disaggregator
```

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

# Consumption
spatial.disagg_households_power(
    by="households",
    weight_by_income=True,
    year=2022,
    scale_by_pop=True,
).to_csv(f"dr_hh_power_demand_2022.csv")

# Timeseries
temporal.disagg_temporal_power_housholds_slp(
    use_nuts3code=True,
    by="households",
    weight_by_income=True,
    year=2022,
    scale_by_pop=True,
).to_csv(f"dr_hh_power_timeseries_2022.csv")
```

### Haushalte: Gas

Zeitreihen je NUTS3

```python
from disaggregator import temporal

# Timeseries
temporal.disagg_temporal_gas_households(
    use_nuts3code=True,
    how='top-down',
    year=2011,
).to_csv(f"dr_hh_gas_timeseries_2011.csv")
```


### GHD und Industrie: Strom

Bedarfe und Zeitreihen je NUTS3
- Bedarfe: Je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: Für alle WZ bedarfsgewichtet aggregiert, Einzelprofile basieren
  je nach WZ auf gemessenen oder SLP inkl. Wanderung
- Letzte verfügbare Daten aus 2017, Fortschreibung für 2022 mit
  Berücksichtigung Beschäftigte und Effizienzgewinne

```python
from disaggregator import spatial, temporal

#######
# CTS #
#######

# Consumption
spatial.disagg_CTS_industry(
    sector='CTS',
    source='power',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_cts_power_demand_2022.csv")
# Timeseries
temporal.disagg_temporal_power_CTS(
    detailed=False,
    use_nuts3code=True,
    year=2022,
).to_csv("dr_cts_power_timeseries_2022.csv")

############
# Industry #
############

# Consumption
spatial.disagg_CTS_industry(
    sector='industry',
    source='power',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_ind_power_demand_2022.csv")
# Timeseries
temporal.disagg_temporal_industry(
    source="power",
    detailed=False,
    use_nuts3code=True,
    no_self_gen=False,
    year=2022,
).to_csv("dr_ind_power_timeseries_2022.csv")
```

### GHD: Gas

Zeitreihen je NUTS3 für alle WZ bedarfsgewichtet aggregiert, Einzelprofile
basieren je nach WZ auf gemessenen oder SLP inkl. Wanderung. Letzte verfügbare
Daten aus 2017, Fortschreibung für 2022 mit Berücksichtigung Beschäftigte und
Effizienzgewinne.

```python
from disaggregator import spatial, temporal

# Timeseries
x=temporal.disagg_temporal_gas_CTS(
    detailed=False,
    use_nuts3code=True,
    year=2011,
).to_csv("dr_cts_gas_timeseries_2011.csv")
```

### Industrie: Gas

Bedarfe und Zeitreihen je NUTS3
- Bedarfe: Je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: Für alle WZ bedarfsgewichtet aggregiert, Einzelprofile basieren
  je nach WZ auf gemessenen oder SLP inkl. Wanderung
- Letzte verfügbare Daten aus 2017, Fortschreibung für 2022 mit
  Berücksichtigung Beschäftigte und Effizienzgewinne

```python
from disaggregator import spatial, temporal

# Consumption
spatial.disagg_CTS_industry(
    sector='industry',
    source='gas',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_ind_gas_demand_2022.csv")
# Timeseries
x=temporal.disagg_temporal_industry(
    source="gas",
    detailed=False,
    use_nuts3code=True,
    no_self_gen=False,
    year=2011,
).to_csv("dr_ind_gas_timeseries_2011.csv")
```
