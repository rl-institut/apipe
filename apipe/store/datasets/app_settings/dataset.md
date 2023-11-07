# Settings für App

Einstellungen für die App.

## Layerliste (rechtes Panel)

- Konfiguration: [config.yml](https://github.com/rl-institut/apipe/blob/main/apipe/store/datasets/app_settings/config.yml) --> `map_panel_layer_list`
- Ergebnisfile: `map_panel_layer_list.json`
- Wird manuell in die App eingepflegt (s.
  [map_config.py](https://github.com/rl-institut/digiplan/blob/main/digiplan/map/map_config.py))

## Settings panels

Die im linken Panel aufgeführten Einstellelemente (Slider und Schalter) werden
hier parametriert.

- Konfiguration des Templates:
  [config.yml](https://github.com/rl-institut/apipe/blob/main/apipe/store/datasets/app_settings/config.yml) --> `panel_settings_templates`
- Parametrierung der Slider und Schalter:
  [panels.py](https://github.com/rl-institut/apipe/blob/main/apipe/store/datasets/app_settings/scripts/panels.py)
- Ergebnisfiles:
    - `energy_settings_panel.json`
    - `heat_settings_panel.json`
    - `traffic_settings_panel.json`
- Werden in die App eingelesen

### Parametrierung der Einstellelemente

Für die Slider werden folgende Attribute gesetzt:
Minimum, Maximum, Schrittweite, Startwert, Status-quo-Wert, Zielwert 2045.
Diese werden wie folgt bestimmt (vgl. auch (i)-Tooltips an den Elementen):

| **Technologie**         | **Element id** | **Maximum**                         | **Startwert**       | **Status-quo-Wert**                             | **Zielwert 2045**                              |
|-------------------------|----------------|-------------------------------------|---------------------|-------------------------------------------------|------------------------------------------------|
| Windenergie             | `s_w_1`        | Inst. Leistung in bestehenden VR/EG | Wie Status-quo-Wert | Inst. Leistung 2022                             | Aus Flächenziel Sachsen-Anhalt (2,2 % in 2032) |
|                         | `s_w_3`        | -                                   | Wie Status-quo-Wert | On                                              | -                                              |
|                         | `s_w_4`        | -                                   | Wie Status-quo-Wert | Off                                             | -                                              |
|                         | `s_w_4_1`      | -                                   | Wie Status-quo-Wert | On                                              | -                                              |
|                         | `s_w_4_2`      | -                                   | Wie Status-quo-Wert | Off                                             | -                                              |
|                         | `s_w_5`        | -                                   | Wie Status-quo-Wert | Off                                             | -                                              |
|                         | `s_w_5_1`      | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
|                         | `s_w_5_2`      | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
| Freiflächen-PV          | `s_pv_ff_1`    |                                     | Wie Status-quo-Wert | Inst. Leistung 2022                             | Aus EEG 2023 und regionalen Potenzialen        |
|                         | `s_pv_ff_3`    | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
|                         | `s_pv_ff_4`    | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
| Aufdach-PV              | `s_pv_d_1`     |                                     | Wie Status-quo-Wert | Inst. Leistung 2022                             | Aus EEG 2023 und regionalen Potenzialen        |
|                         | `s_pv_d_3`     | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
|                         | `s_pv_d_4`     | 100 %                               | Wie Status-quo-Wert | Aus MaStR                                       | -                                              |
| Wasserkraft             | `s_h_1`        | Inst. Leistung 2022                 | Wie Status-quo-Wert | Inst. Leistung 2022                             | Inst. Leistung 2022                            |
| Stromverbrauch          | `s_v_1`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `s_v_3`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `s_v_4`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `s_v_5`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
| Batterie-Großspeicher   | `s_s_g_1`      | 50 %                                | Wie Status-quo-Wert | Aus inst. Kapazität und Einspeisung 2022        | -                                              |
|                         |                |                                     | Wie Status-quo-Wert |                                                 |                                                |
| WP dezentral            | `w_d_wp_1`     | 95 %                                | 50 %                | Inst. Leistung 2022 aus BMWK Langfristszenarien | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_d_wp_3`     | 95 %                                | 50 %                | -                                               | -                                              |
|                         | `w_d_wp_4`     | 95 %                                | 50 %                | -                                               | -                                              |
|                         | `w_d_wp_5`     | 95 %                                | 50 %                | -                                               | -                                              |
| WP zentral              | `w_z_wp_1`     | 95 %                                | 50 %                | Inst. Leistung 2022 aus BMWK Langfristszenarien | Wert 2045 aus BMWK Langfristszenarien          |
| Wärmeverbrauch          | `w_v_1`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_v_3`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_v_4`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_v_5`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
| Wärmespeicher dezentral | `w_d_s_1`      | 200 %                               | 100 %               | -                                               | -                                              |
| Wärmespeicher zentral   | `w_z_s_1`      | 200 %                               | 100 %               | -                                               | -                                              |

Die Maxima der Regler im Hauptpanel (`s_w_1`, `s_pv_ff_1` usw.) werden in der
App dynamisch aus den durch die UserInnen vorgenommenen Detaileinstellungen
(`s_w_3`, `s_pv_ff_1` usw.) berechnet.
