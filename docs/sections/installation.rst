.. _installation_label:

Installation
============

**Hinweis: Nur Linux, Windows wird derzeit nicht unterstützt.**

Klonen Sie zunächst über SSH mit

.. code-block::

   git clone git@github.com:rl-institut-private/digipipe.git /local/path/to/digipipe/


Installation via pip
--------------------

Stellen Sie sicher, dass Sie Python >= 3.6 installiert haben.
Anschließend erstellen Sie eine virtuelle Umgebung mit:

.. code-block::

   virtualenv --python=python3.8 venv
   source venv/bin/activate


Es werden einige zusätzliche Systempakete benötigt, die Sie mit dem folgenden Befehl installieren

.. code-block::

   apt install gdal-bin python3-gdal libspatialindex-dev imagemagick osmium-tool graphviz graphviz-dev

Weitere Anmerkungen:

* Stellen Sie sicher, dass Sie `GDAL>=3.0` haben, da ältere Versionen nicht
  funktionieren.
* `imagemagick` ist optional und nur für die Berichtserstellung erforderlich

Installieren Sie das Paket mit:

.. code-block::

   pip install -e /local/path/to/digipipe/


Installation via conda
----------------------

Stellen Sie sicher, dass Sie conda installiert haben
(bspw. miniconda). Erstellen Sie dann die env:

.. code-block::

   conda create -n digipipe /local/path/to/digipipe/environment.yml
   conda activate digipipe



