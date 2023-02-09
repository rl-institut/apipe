import copy
import os

import json
import pandas as pd


here = os.path.dirname(__file__)


class B3Schema:
    def __init__(self, index, columns):
        self.index = index
        self.columns = columns

    @classmethod
    def load_from_csv(cls, path):
        df = pd.read_csv(path, delimiter=";")

        df.columns.name = "field"

        df.index = ["type", "description"]

        index = df.iloc[:, 0]

        columns = df.iloc[:, 1:]

        return cls(index, columns)


SCHEMA_SCAL = B3Schema.load_from_csv(os.path.join(here, "scalars.csv"))

SCHEMA_TS = B3Schema.load_from_csv(os.path.join(here, "timeseries.csv"))
