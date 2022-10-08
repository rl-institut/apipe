"""
Snakefile for this dataset

The file will be automatically detected and included in the Snakemake workflow.
"""

#from zipfile import ZipFile

from digipipe.store.utils import get_abs_dataset_path

#configfile: get_abs_dataset_path("1_preprocessed", "osm_filtered", data_dir=False) / "config.yml"

rule bkg_vg250_extract:
    input:
        get_abs_dataset_path("0_raw", "bkg_vg250") /
        "vg250_01-01.utm32s.gpkg.ebenen.zip"
    output:
        get_abs_dataset_path("1_preprocessed", "bkg_vg250") / "bkg_vg250.gpkg"
    params:
        outpath=get_abs_dataset_path("1_preprocessed", "bkg_vg250"),
        original_file=get_abs_dataset_path("1_preprocessed", "bkg_vg250") / "DE_VG250.gpkg",
        file_path_in_zip=str("vg250_01-01.utm32s.gpkg.ebenen/vg250_ebenen_0101/DE_VG250.gpkg")
    shell:
        """
        unzip -j {input} {params.file_path_in_zip} -d {params.outpath}
        ogr2ogr -f GPKG -t_srs EPSG:3035 {output} {params.original_file}
        rm {params.original_file}
        """
        # with ZipFile(input) as f:
        #     f.extract(
        #         "vg250_01-01.utm32s.gpkg.ebenen/vg250_ebenen_0101/DE_VG250.gpkg",
        #         path=get_abs_dataset_path("1_preprocessed", "bkg_vg250")
        #     )
        #     #f.extractall(get_abs_dataset_path("1_preprocessed", "bkg_vg250"))
        #     #f.extractall()
