import os
import json


def generate_dataset_mds():
    """Create markdown files for each dataset category

    This function generates/updates Markdown files based on the content of
    'dataset.md' and 'metadata.json' in each dataset folder per category.
    One Markdown file per category is generated (containing data for every
    dataset of this category).
    The Markdown files are saved to the 'docs/datasets' directory.

    Returns
    -------
    None
    """
    source_dir = "digipipe/store"
    target_dir = "docs/datasets"

    categories = ["raw", "preprocessed", "datasets"]

    for category in categories:
        category_dir = os.path.join(source_dir, category)
        category_md_file = os.path.join(target_dir, f"{category}_datasets.md")

        with open(category_md_file, "w") as md_file:
            # Write category heading
            md_file.write(f"# '{category.capitalize()}' Datasets \n")

            # Traverse the category directory
            for root, dirs, files in os.walk(category_dir):
                if "TEMPLATE" in dirs:
                    dirs.remove("TEMPLATE")
                for file in files:
                    if file == "dataset.md":
                        dataset_md_file = os.path.join(root, file)
                        # Write the content of the dataset.md file
                        with open(dataset_md_file, "r") as dataset_file:
                            md_file.write("\n------------------------------\n")
                            md_file.write("#" + dataset_file.read())
                            md_file.write("\n")
                            md_file.write(
                                f"**Dataset: `{category}/{os.path.basename(root)}`**\n\n"
                            )

                        # Check if the corresponding metadata.json file exists
                        metadata_json_file = os.path.join(root, "metadata.json")
                        if os.path.exists(metadata_json_file):
                            # Read the metadata from metadata.json
                            with open(metadata_json_file, "r") as metadata_file:
                                metadata = json.load(metadata_file)
                            # Write the metadata section
                            md_file.write('??? metadata "Metadata"\n')
                            md_file.write("    ```json\n")
                            md_file.write(
                                "    "
                                + json.dumps(metadata, indent=4).replace(
                                    "\n", "\n    "
                                )
                            )
                            md_file.write("\n    ```\n")

        print(f"Generated {category}.md file.")

    print("Generation of dataset category markdown files completed.")


generate_dataset_mds()
