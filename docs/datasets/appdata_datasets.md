# 'Appdata' Datasets 

------------------------------
## Datapackage für App

Von der App benötigte Daten in JSON `datapackage_app.json`.

Generelle Struktur
- `<KATEGORIE>`
  - `<DATENTYP>` (`scalars`, `sequences` oder `geodata`)
    - `<DATENSATZNAME>`
      - `description`: Beschreibung
      - `path`: Pfad zur Datei im Datapackage
      - `fields`: Felder-/Spaltendefinition
      - `_source_path`: Pfad zur Datei im Quelldatensatz
        - `dataset`
        - `file`

Kategorien bzw. Inhalt `resources`:
- `base_data`: Basisdaten
- `production`: Energiebereitstellung
- `demand`: Energiebedarf
- `potential_areas`: Potenzialgebiete EE

**Dataset: `appdata/datapackage`**


------------------------------
## Settings für App

Einstellungen für die App.

### Layerliste (rechtes Panel)

- Konfiguration: [config.yml](config.yml) -> `map_panel_layer_list`
- Ergebnisfile: `map_panel_layer_list.json`
- Wird manuell in die App eingepflegt (s. [map_config.py](https://github.com/rl-institut-private/digiplan/blob/dev/digiplan/map/map_config.py))

### Settings panels

- Konfiguration des Templates: [config.yml](config.yml) -> `panel_settings_templates`
- Ergebnisfiles:
  - `energy_settings_panel.json`
  - `heat_settings_panel.json`
  - `traffic_settings_panel.json`
- Werden in die App eingelesen

**TODO**: Parametrierung der Slider & Switches beschreiben

**Dataset: `appdata/settings`**

