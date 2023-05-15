import os
import yaml

# --- WORK IN PROGRESS ---
#
# Define the categories and their corresponding directories
categories = {
    "Raw": "raw",
    "Preprocessed": "preprocessed",
    "Datasets": "datasets",
    "Appdata": "appdata",
}

# Traverse the directory tree and collect the datasets
datasets = {}
for category, directory in categories.items():
    datasets[category] = []
    for root, dirs, files in os.walk(f"digipipe/store/{directory}"):
        for file in files:
            if file.endswith(".mds"):
                path = os.path.join(root, file)
                with open(path, "r") as f:
                    content = f.read()
                    datasets[category].append(
                        {"name": file[:-4], "description": content}
                    )

# Modify the mkdocs.yml file
with open("mkdocs.yml", "r") as f:
    config = yaml.safe_load(f)

nav = config.get("nav", [])
for category, datasets in datasets.items():
    nav.append(
        {
            category: [
                f"{dataset['name']}: digipipe/store/{category}/{dataset['name']}.md"
                for dataset in datasets
            ]
        }
    )

config["nav"] = nav

with open("mkdocs.yml", "w") as f:
    yaml.dump(config, f)
