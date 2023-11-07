# Erzeugungsanlagen aus Marktstammdatenregister

Ereugungsanlagen aus dem Markstammdatenregister, das mit dem Tool
[open-mastr](https://github.com/OpenEnergyPlatform/open-MaStR) erstellt und
abgelegt wurde. Die Daten wurden folgendermaßen erstellt:
```
from open_mastr import Mastr
db = Mastr()
db.download("bulk")
db.to_csv(None)  # (None for all data)
```

Die abgelegten CSV-Dateien (alle Tabellen) wurden um einen benutzerdefinierten
Export von Speichereinheiten mit
`sqlite3 -header -csv -separator "," open-mastr.db "select * from storage_units;" > bnetza_mastr_storage_unit_raw.csv`
erweitert. Anschließend wurden alle Dateien komprimiert.

Das Marktstammdatenregister (MaStR) ist ein deutsches Register, welches von der
Bundesnetzagentur (BNetza) bereitgestellt wird und alle in Deutschland
befindlichen Strom- und Gasanlagen erfasst.
