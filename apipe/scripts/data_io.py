import json
import os
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests


def load_json(file_path: Path) -> dict:
    with open(file_path, "r") as f:
        return json.load(f)


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
    # Loop through each directory and copy files within `data` directory
    for d in dir_list:
        src_dir = os.path.join(src_path, d, "data")
        dst_dir = os.path.join(dest_path, d, "data")
        # Check if destination directory for data already exists
        if os.path.exists(dst_dir):
            files_to_copy = []
            for file in os.listdir(src_dir):
                src_file = os.path.join(src_dir, file)
                dst_file = os.path.join(dst_dir, file)
                if os.path.isfile(dst_file):
                    if file == ".gitkeep":
                        continue
                    print(f"\n'{file}' already exists in '{d}/data'.")
                    overwrite_file = input("Do you want to update it? (y/n) ")
                    while overwrite_file.lower() not in ["y", "n"]:
                        overwrite_file = input(
                            """Invalid input. Enter 'y' or 'n': """
                        )
                    if overwrite_file.lower() == "y":
                        shutil.copy(src_file, dst_file)
                        print(f"'{file}' updated.")
                    else:
                        continue
                else:
                    files_to_copy.append(file)

            for file in files_to_copy:
                src_file = os.path.join(src_dir, file)
                dst_file = os.path.join(dst_dir, file)
                shutil.copy(src_file, dst_file)
                print(f"\n'{file}' added to '{d}/data'.")
        # If directory doesn't exist, copy it entirely
        else:
            shutil.copytree(src_dir, dst_dir)
            print(f"\nAdded '{d}' to 'store/raw'.")


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


def df_reduce_mem_usage(
    df: pd.DataFrame, show_reduction: bool = False
) -> pd.DataFrame:
    """Function to automatically check if columns of a pandas DataFrame can
    be reduced to a smaller data type. Source:
    https://www.mikulskibartosz.name/how-to-reduce-memory-usage-in-pandas/

    Parameters
    ----------
    df: pd.DataFrame
        DataFrame to reduce memory usage on
    show_reduction : bool
        If True, print amount of memory reduced

    Returns
    -------
    pd.DataFrame
        DataFrame with memory usage decreased
    """
    start_mem = df.memory_usage().sum() / 1024 ** 2

    for col in df.columns:
        col_type = df[col].dtype

        if col_type != object and str(col_type) != "category":
            c_min = df[col].min()
            c_max = df[col].max()

            if str(col_type)[:3] == "int":
                if (
                    c_min > np.iinfo(np.int16).min
                    and c_max < np.iinfo(np.int16).max
                ):
                    df[col] = df[col].astype("int16")
                elif (
                    c_min > np.iinfo(np.int32).min
                    and c_max < np.iinfo(np.int32).max
                ):
                    df[col] = df[col].astype("int32")
                else:
                    df[col] = df[col].astype("int64")
            else:
                if (
                    c_min > np.finfo(np.float32).min
                    and c_max < np.finfo(np.float32).max
                ):
                    df[col] = df[col].astype("float32")
                else:
                    df[col] = df[col].astype("float64")

        else:
            df[col] = df[col].astype("category")

    end_mem = df.memory_usage().sum() / 1024 ** 2

    if show_reduction is True:
        print(
            "Reduced memory usage of DataFrame by "
            f"{(1 - end_mem/start_mem) * 100:.2f} %."
        )

    return df
