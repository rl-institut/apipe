import os
import shutil
import zipfile

import requests


def download_file(url: str, dest_file: str) -> None:
    """Downloads file from url to specified destination.

    Parameters
    ----------
    url : str
        The URL string of the file to be downloaded.
    dest_file : str
        The destination path/filename of the file to be downloaded.
    """
    print(f"Downloading from {url}")
    with requests.get(url, stream=True, timeout=10) as r:
        r.raise_for_status()
        with open(dest_file, "wb") as f:
            for i, block in enumerate(r.iter_content(chunk_size=1024)):
                f.write(block)
                downloaded_size = i * 1024 / 1e6
                print(f"{downloaded_size:.2f} MB...", end="\r")
    print(f"Download completed to {dest_file}")


def extract_zipfile(zip_path: str, dest_path: str) -> None:
    """Extracts zipfile to destination path

    Parameters
    ----------
    zip_path : str
        Path to the zipfile to be extracted
    dest_path : str
        Destination path where the contents of zipfile are to be extracted
    """
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(dest_path)


def copy_files(src_path: str, dest_path: str) -> None:
    """Copies files from source to destination directory

    Parameters
    ----------
    src_path : str
        Source path of directory to copy from
    dest_path : str
        Destination path of directory to copy to
    """
    # Get a list of directories within source directory
    dir_list = [
        d
        for d in os.listdir(src_path)
        if os.path.isdir(os.path.join(src_path, d))
    ]
    # Loop through each directory and copy files within "data" directory to
    # destination directory of dataset
    for d in dir_list:
        src_dir = os.path.join(src_path, d, "data")
        dst_dir = os.path.join(dest_path, d, "data")
        # Check if destination directory for dataset already exists
        if os.path.exists(dst_dir):
            print(f"\nRaw data for '{d}' already exists. ")
            overwrite = input("Do you want to update it? (y/n) ")
            while overwrite.lower() not in ["y", "n"]:
                overwrite = input("Invalid input. Enter 'y' or 'n': ")
            if overwrite.lower() == "y":
                for file in os.listdir(src_dir):
                    src_file = os.path.join(src_dir, file)
                    shutil.copy(src_file, dst_dir)
                print(f"Updated data for '{d}'.")
            else:
                continue
        else:
            shutil.copytree(src_dir, dst_dir)
            print(f"\nAdded new dataset '{d}'.")


def clean_folder(folder_path: str) -> None:
    """Cleans up a given folder by removing all files and directories within it
    (except for .gitkeep)

    Parameters
    ----------
    folder_path : str
        Path to the folder that should be cleaned
    """
    # Clean temp folder
    folder_content = os.listdir(folder_path)
    for item in folder_content:
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path) and item != ".gitkeep":
            os.remove(item_path)
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)
