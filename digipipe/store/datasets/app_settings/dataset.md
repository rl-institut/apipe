# Settings für App

Einstellungen für die App.

## Layerliste (rechtes Panel)

- Konfiguration: [config.yml](config.yml) --> `map_panel_layer_list`
- Ergebnisfile: `map_panel_layer_list.json`
- Wird manuell in die App eingepflegt (s. [map_config.py](https://github.com/rl-institut-private/digiplan/blob/dev/digiplan/map/map_config.py))

## Settings panels

- Konfiguration des Templates: [config.yml](config.yml) --> `panel_settings_templates`
- Ergebnisfiles:
    - `energy_settings_panel.json`
    - `heat_settings_panel.json`
    - `traffic_settings_panel.json`
- Werden in die App eingelesen

**TODO**: Parametrierung der Slider & Switches beschreiben

- `s_pv_d_1`: Installierbare Leistung PV-Aufdachanlagen.
  Max. 50 % aller Dächer von nicht-denkmalgeschützten Gebäuden mit Ausrichtung
  Süden, Osten, Westen und Flachdächern.
