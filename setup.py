from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="digipipe",
    version="0.0.0",
    description="Pipeline for data and energy system in the Digiplan project.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rl-institut/digipipe/",
    author="Reiner Lemoine Institut",
    author_email='jonathan.amme@rl-institut.de',
    license='GNU AGPLv3',
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        # Indicate who your project is intended for
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        # Pick your license as you wish
        "License :: OSI Approved :: AGPL License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    packages=find_packages(),
    python_requires=">=3, <4",
    install_requires=requirements,
    #extras_require={"dev": [], "test": []},  # Optional
    project_urls={  # Optional
        "Bug Reports": "https://github.com/rl-institut/digipipe/issues",
        "Source": "https://github.com/rl-institut/digipipe",
    },
)
