import json
from pathlib import Path

# Gegebene Listen
region_list = [
    "r120640428428",
    "r120640472472",
    "r120670124124",
    "r120670201201",
]
file_list = [
    "residual_waste-bpchp_heat_low_central",
    "electricity-heatpump_decentral",
    "biomass_solid-bpchp_heat_low_decentral",
    "biomass_gas-bpchp_heat_low_decentral",
    "electricity-heatpump_central",
    "electricity-electrolyzer",
    "electricity-transmission",
    "heat_low_decentral-solarthermal_plant",
    "heat_low_central-storage",
    "biomass_solid-commodity",
    "residual_waste-commodity",
    "biomass_gas-bpchp_heat_high",
    "electricity-wind",
    "h2-bpchp_heat_high",
    "electricity-boiler_heat_high",
    "heat_low_decentral-demand_ind",
    "heat_low_decentral-demand_hh",
    "h2-export",
    "heat_low_central-demand_hh",
    "heat_low_central-demand_ind",
    "electricity-pv_rooftop",
    "bus",
    "heat_low_decentral-storage",
    "heat_low_central-demand_cts",
    "h2-commodity",
    "electricity-pv_agri_horizontal",
    "heat_low_decentral-demand_cts",
    "biomass_gas-bpchp_heat_low_central",
    "electricity-small_battery_storage",
    "biomass_gas-commodity",
    "electricity-demand_mob",
    "electricity-heatpump_heat_high",
    "biomass_solid-bpchp_heat_high",
    "electricity-demand_ind",
    "biomass_solid-bpchp_heat_low_central",
    "electricity-pv_agri_vertical",
    "residual_waste-bpchp_heat_high",
    "biomass_solid-boiler_heat_high",
    "electricity-large_battery_storage",
    "electricity-export",
    "h2-boiler_heat_high",
    "heat_high-demand_ind",
    "electricity-import",
    "electricity-pv_ground",
    "electricity-demand_hh",
    "electricity-demand_cts",
    "h2-transmission",
]

# Filter für CapacityCost (ohne "commodity", "demand", "bus", "export", "import", "transmission")
capacity_cost_files = [
    f
    for f in file_list
    if not any(
        x in f
        for x in [
            "commodity",
            "demand",
            "bus",
            "export",
            "import",
            "transmission",
        ]
    )
]

# Filter für Demand ("demand" enthalten)
demand_files = [f for f in file_list if "demand" in f]

# Filter für MarginalCost ("commodity", "transmission", "import", "export" enthalten)
marginal_cost_files = [
    f
    for f in file_list
    if any(x in f for x in ["commodity", "transmission", "import", "export"])
]


cost_perturbations = {
    f"Perturbation{i + 1}": [
        {
            "FileName": filename,
            "Regions": region_list,
            "Column": "capacity_cost",
            "PerturbationMethod": "multiplication",
            "PerturbationParameter": [0.5, 0.9, 1.2, 1.5],
        }
    ]
    for i, filename in enumerate(capacity_cost_files)
}

# Marginal Cost als separate Perturbationen
for i, filename in enumerate(
    marginal_cost_files, start=len(capacity_cost_files) + 1
):
    cost_perturbations[f"Perturbation{i}"] = [
        {
            "FileName": filename,
            "Regions": region_list,
            "Column": "marginal_cost",
            "PerturbationMethod": "addition",
            "PerturbationParameter": [],  # Werte manuell einfügen
        }
    ]

# Demand-Perturbationen als eigene Perturbations
demand_perturbations = {
    f"Perturbation{i + 1}": [
        {
            "FileName": filename,
            "Regions": region_list,
            "Column": "amount",
            "PerturbationMethod": "multiplication",
            "PerturbationParameter": [0.5, 1.2, 1.5, 2.0],
        }
    ]
    for i, filename in enumerate(demand_files)
}

# Gesamtstruktur des JSON
json_data = {
    "Preferences": {
        "SettingName": "ScenarioAnalysis",
        "NumberOfPerturbations": 2,
        "BaseScenName": "2045_scenario",
    },
    "CostPerturbations": cost_perturbations,
    "DemandPerturbations": demand_perturbations,
}

# JSON in eine Datei speichern
output_path = Path(Path(__file__).parent, "PerturbationJSON.json")
with open(output_path, "w") as f:
    json.dump(json_data, f, indent=2)
