"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path


rule extract:
    input:
        get_abs_dataset_path("raw", "bkg_vg250") /
        "vg250_01-01.utm32s.gpkg.ebenen.zip"
    output:
        get_abs_dataset_path("preprocessed", "bkg_vg250") / "bkg_vg250.gpkg"
    params:
        outpath=get_abs_dataset_path("preprocessed", "bkg_vg250"),
        original_file=get_abs_dataset_path("preprocessed", "bkg_vg250") / "DE_VG250.gpkg",
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
        #         path=get_abs_dataset_path("preprocessed", "bkg_vg250")
        #     )
        #     #f.extractall(get_abs_dataset_path("preprocessed", "bkg_vg250"))
        #     #f.extractall()
