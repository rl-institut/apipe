import os
import shutil
import json


def copy_dataset_md_files():
    source_dir = "digipipe/store"
    target_dir = "docs/datasets"

    sub_dirs = ["raw", "preprocessed", "datasets"]

    for sub_dir in sub_dirs:
        root_dir = os.path.join(source_dir, sub_dir)
        target_sub_dir = os.path.join(target_dir, sub_dir)
        os.makedirs(target_sub_dir, exist_ok=True)

        for root, dirs, files in os.walk(root_dir):
            # Extrahieren des entsprechenden Sub-Verzeichnisnamens
            sub_dir_path = root.split(root_dir + os.sep, 1)
            if len(sub_dir_path) == 2:
                sub_dir_path = sub_dir_path[1]
            else:
                continue

            # Durchsuchen der "dataset.md" und "metadata.json"-Dateien im aktuellen Verzeichnis
            dataset_md_file = None
            metadata_json_file = None
            for file in files:
                if file == "dataset.md":
                    dataset_md_file = os.path.join(root, file)
                elif file == "metadata.json":
                    metadata_json_file = os.path.join(root, file)

            # Wenn eine "dataset.md"-Datei gefunden wurde
            if dataset_md_file and sub_dir_path != ".TEMPLATE":
                # Kopieren der "dataset.md"-Datei in das Zielverzeichnis
                target_dataset_md_file = os.path.join(
                    target_sub_dir, sub_dir_path + ".md"
                )
                if os.path.isfile(
                    dataset_md_file
                ):  # Überprüfung, ob die Datei existiert
                    shutil.copy2(dataset_md_file, target_dataset_md_file)

                # Wenn eine "metadata.json"-Datei vorhanden ist, den Inhalt als Markdown-Codeblock einfügen
                if metadata_json_file:
                    # Laden des Inhalts der "metadata.json"-Datei
                    with open(metadata_json_file, "r") as json_file:
                        metadata = json.load(json_file)

                    # Öffnen der "datasets.md"-Datei im Zielverzeichnis im Anhang-Modus
                    with open(target_dataset_md_file, "a") as dataset_md:
                        # Hinzufügen des JSON-Inhalts als Markdown-Codeblock am Ende der Datei
                        dataset_md.write("\n\n")
                        dataset_md.write('??? info "Metadaten"\n')
                        dataset_md.write("    ```json\n")
                        dataset_md.write(
                            "    "
                            + json.dumps(metadata, indent=4).replace(
                                "\n", "\n    "
                            )
                        )
                        dataset_md.write("\n    ```\n")

    print(
        "Kopieren der dataset.md-Dateien und Einfügen von metadata.json abgeschlossen."
    )


copy_dataset_md_files()
