# 'Appdata' Datasets 

------------------------------
## Datapackage für App

Von der App benötigte Daten in JSON `datapackage_app.json`.

Generelle Struktur:

- `<KATEGORIE>`
  - `<DATENTYP>` (`scalars`, `sequences` oder `geodata`)
    - `<DATENSATZNAME>`
      - `description`: Beschreibung
      - `path`: Pfad zur Zieldatei im Datapackage
      - `fields`: Felder-/Spaltendefinition
      - `_source_path`: Pfad zur Datei im Quelldatensatz
        - `dataset`: Name in `store/datasets/`
        - `file`: Datei

Kategorien bzw. Inhalt `resources`:

- `base_data`: Basisdaten
- `production`: Energiebereitstellung
- `demand`: Energiebedarf
- `potential_areas`: Potenzialgebiete EE
- `emissions`: Emissionen
- `settings`: App-Settings
- `captions`: App-Captions

**Dataset: `appdata/datapackage`**

