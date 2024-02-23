# Potenzialgebiete Agri-PV (Region)

## Methodik

### Add dataset for region: "potentialarea_pv_ground_region2"

- Clip to region
- Vectorize
- Delete areas below min. area (in ha (1 ha = 1 cell)) threshold (from config, default setting=0)
- Calculate average pixel value and assign to polygon

### Result files (geodata)

- potentialarea_pv_ground_soil_quality_low_region.gpkg (1)
- potentialarea_pv_ground_soil_quality_medium_region.gpkg (2)
- potentialarea_pv_ground_permanent_crops_region.gpkg (3)

### Result files (statistics, see dataset "potentialarea_pv_ground_region" in digipipe)

- per mun:
potentialarea_pv_ground_area_stats_muns.csv
- total shares:
potentialarea_pv_ground_area_shares.json
- targets: potentialarea_pv_ground_regionalized_targets.json
