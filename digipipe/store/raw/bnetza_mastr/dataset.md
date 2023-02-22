# Raw dataset: bnetza_mastr

Power units from Marktstammdatenregister obtained and dumped with the
[open-mastr](https://github.com/OpenEnergyPlatform/open-MaStR) tool. The
dump was created with:

```
from open_mastr import Mastr
db = Mastr()
db.download("bulk")
db.to_csv(None)  # (None for all data)
```

The dumped CSV files (all tables) have been extended by a custom export of
storage units with
`sqlite3 -header -csv -separator "," open-mastr.db "select * from storage_units;" > bnetza_mastr_storage_unit_raw.csv`.
Subsequently, all files were zipped.

The Marktstammdatenregister (MaStR) is a German register provided by the
German Federal Network Agency (Bundesnetzagentur / BNetza) that keeps
track of all power and gas units located in Germany.
