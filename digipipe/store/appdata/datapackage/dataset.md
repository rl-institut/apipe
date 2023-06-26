# Datapackage für App

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
