import builtins
import os
import shutil
import zipfile


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
    # Get a list of all directories within source directory
    dir_list = [
        d
        for d in os.listdir(src_path)
        if os.path.isdir(os.path.join(src_path, d))
    ]
    # Loop through each directory and copy it to destination directory
    for d in dir_list:
        src_dir = os.path.join(src_path, d)
        dst_dir = os.path.join(dest_path, d)
        if os.path.exists(dst_dir):
            print(
                f"Directory '{dst_dir}' already exists."
                "Do you want to overwrite it?"
            )
            overwrite = builtins.input("Enter y/n: ")
            if overwrite.lower() == "y":
                shutil.rmtree(dst_dir)
            else:
                continue
        shutil.copytree(src_dir, dst_dir)


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
