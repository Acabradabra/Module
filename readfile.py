# -*- coding: utf-8 -*-
from pathlib import Path
from numpy.typing import NDArray
import mmap
import numpy as np

from time import perf_counter
t0 = perf_counter()

class FluentInterpolationParser:
    def __init__(self, fname: str | Path, float_type: type = float):
        if not Path(fname).exists():
            raise FileNotFoundError(f"Fluent IP file not found {fname}")

        self._names = []

        with open(fname, "rb") as fp:
            self._parse_manager(float_type, fp)

        # TODO maybe working on 1-D array and reshaping would be
        # faster `data = data.reshape((num_points, num_cols))`?
        self._data = self._data.T

    def _parse_manager(self, float_type, fp):
        """ Handle memory mapping with read-only access for parsing. """
        with mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            parser = iter(mm.readline, b"")

            self._n_first      = int(next(parser).decode())
            self._n_dimensions = int(next(parser).decode())
            self._n_points     = int(next(parser).decode())
            self._n_variables  = int(next(parser).decode())

            for _ in range(self._n_variables):
                self._names.append(next(parser).decode().strip("\n"))

            num_cols = self._n_variables + self._n_dimensions
            shape = (num_cols, self._n_points)

            self._data = np.empty(shape, float_type)
            self._parse_data(float_type, parser)

    def _parse_data(self, float_type, parser):
        """ Manage parsing of all variables from raw text. """
        vars_finish = self._n_variables + 5

        for var_idx in range(self._n_variables):
            start_idx = vars_finish + var_idx * (self._n_points + 1)
            end_idx   = start_idx + self._n_points
            print(f"var = {var_idx:>3} ({start_idx:>10}:{end_idx:>10})")
            self._parse_variable(float_type, parser, self._data[var_idx])

    def _parse_variable(self, float_type, parser, data):
        """ Manage parsing of a single variable into an array segment. """
        for pnt_idx in range(self._n_points):
            data[pnt_idx] = float_type(next(parser).decode().strip("("))

        # XXX: ignore `)` in the end so that parser is ready for next:
        _ = next(parser) 

    @property
    def n__dimensions(self) -> int:
        """ Provides access to number of dimensions. """
        return self._n_dimensions
    
    @property
    def n_points(self) -> int:
        """ Provides access to number of points. """
        return self._n_points

    @property
    def n_variables(self) -> int:
        """ Provides access to number of variables. """
        return self._n_variables

    @property
    def variable_names(self) -> list[str]:
        """ Provides access to list of variables names. """
        return self._names
    
    @property
    def data(self) -> NDArray[np.float64]:
        """ Return access to the whole data table. """
        return self._data
    
    def get_data(self, name) -> NDArray[np.float64]:
        """ Retrieve data for selected variable name. """
        if name not in self._names:
            raise KeyError(f"Unknown variable `{name}`")
        return self._data[:, self._names.index(name)]


# workflow = FluentInterpolationParser("sample.ip")
# print(f"Took {perf_counter()-t0} s")
