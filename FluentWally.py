# -*- coding: utf-8 -*-
from io import StringIO
from pathlib import Path
from numpy.typing import NDArray
from pyparsing import OneOrMore
from pyparsing import Word
from pyparsing import Suppress
import mmap
import numpy as np
import pandas as pd
import pyparsing as pp


class FluentInterpolationParser:
    """ Provides parsing of exported Fluent interpolation files.

    Parameters
    ----------
    fname: str | Path
        Path to `.ip` file to be parsed; file must be stored in plain
        text format, by default in UTF-8 encoding.
    float_type: type = float
        Floatting point type; using built-in is faster, but you may also
        use some NumPy type at your convenience (much slower).
    encoding: str = "utf-8"
        Encoding type for reading the file.
    """
    def __init__(self, fname: str | Path, float_type: type = float,
                 encoding: str = "utf-8"):
        if not Path(fname).exists():
            raise FileNotFoundError(f"Fluent IP file not found {fname}")

        with open(fname, "r", encoding=encoding) as fp:
            self._parse_manager(float_type, fp)

    def _parse_manager(self, float_type, fp):
        """ Handle memory mapping with read-only access for parsing. """
        with mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            parser = iter(mm.readline, b"")

            self._n_first      = int(next(parser).decode())
            self._n_dimensions = int(next(parser).decode())
            self._n_points     = int(next(parser).decode())
            self._n_variables  = int(next(parser).decode())
            self._n_columns    = self._n_variables + self._n_dimensions

            self._names = [f"x{i}" for i in range(self._n_dimensions)]

            for _ in range(self._n_variables):
                self._names.append(next(parser).decode()\
                                   .strip("\r").strip("\n").strip("\r\n"))

            shape = (self._n_columns, self._n_points)

            self._data = np.empty(shape, float_type)
            self._parse_data(float_type, parser)

        # TODO maybe working on 1-D array and reshaping would be
        # faster `data = data.reshape((num_points, num_cols))`?
        self._data = self._data.T

    def _parse_data(self, float_type, parser):
        """ Manage parsing of all variables from raw text. """
        vars_finish = self._n_variables + 5

        for var_idx in range(self._n_columns):
            name = self._names[var_idx]
            s_idx = vars_finish + var_idx * (self._n_points + 1)
            e_idx = s_idx + self._n_points
            print(f"{name:>15s} in range ({s_idx:>10}:{e_idx:>10})")

            self._parse_variable(float_type, parser, self._data[var_idx])

    def _parse_variable(self, float_type, parser, data):
        """ Manage parsing of a single variable into an array segment. """
        for pnt_idx in range(self._n_points):
            data[pnt_idx] = float_type(next(parser).decode().strip("("))

        # XXX: ignore `)` in the end so that parser is ready for next:
        _ = next(parser)

    @property
    def n_dimensions(self) -> int:
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


class FluentInputRow:
    """ Fluent TSV file data entry for automated parameter generation.

    TODO
    ====
    - Clarify the meaning of `parameterid` in Ansys Fluent.
    - Allow operations between input variables (overloads).

    Parameters
    ==========
    name: str
        Name of parameter, no spaces allowed.
    value: str
        String representation of the value, if required with units.
    unit: str = ""
        Units of parameter.
    kind: str = ""
        Unit family of parameter.
    inpar: bool = False
        Whether to use parameter as a design input parameter.
    outpar: bool = False
        Whether to use parameter as a design output parameter.
    descr: str = ""
        A short description of the parameter.
    """
    def __init__(self,
            name: str,
            value: str,
            unit: str = "",
            kind: str = "",
            inpar: bool = False,
            outpar: bool = False,
            descr: str = ""
        ):
        # Handle values provided with no units:
        if unit and unit not in value:
            value = f"{value} {unit}"

        self._name   = name
        self._value  = value
        self._descr  = descr

        # If a unit family (kind) is provided, use it as unit:
        self._unit   = kind if kind else unit

        # Handle specific flags formatting:
        self._inpar  = self._parflag(inpar)
        self._outpar = self._parflag(outpar)

    def _parflag(self, test):
        return "#t" if test else "#f"

    def __repr__(self):
        return "\t".join((
            f'"{self._name}"',   # name
            f'"{self._value}"',  # definition
            f'"{self._descr}"',  # description
            f'""',               # parameterid
            f'"{self._name}"',   # parametername
            f'"{self._unit}"',   # unit
            f'{self._inpar}',    # input-parameter
            f'{self._outpar}',   # output-parameter
        ))

    @property
    def name(self):
        """ Provides access to entry name. """
        return self._name


class FluentInputFile:
    """ Representation of an Ansys Fluent expressions input file.

    Attributes
    ==========
    HEAD : tuple[str]
        Tuple containing TSV file headers.
    """
    HEAD = ("name", "definition", "description", "parameterid",
            "parametername", "unit", "input-parameter",
            "output-parameter")

    def __init__(self,
            rows: list[FluentInputRow],
            sort_list: bool = False
        ) -> None:
        if sort_list:
            rows = sorted(rows, key=lambda x: x.name)

        self._data = ["\t".join(self.HEAD), *(map(repr, rows))]

    def __str__(self):
        return "\t\n".join((*self._data, ""))

    def append(self, row: FluentInputRow):
        """ Append a single row to the inputs file. """
        self._data.append(repr(row))

    def to_file(self, saveas: str, overwrite: bool = False):
        """ Generate Fluent input file. """
        if Path(saveas).exists() and not overwrite:
            print(f"File {saveas} exists, use `overwrite=True`.")
            return

        with open(saveas, "w") as fp:
            fp.write(str(self))


class FluentSchemePatch:
    """ Parse patch name from a Scheme file. """
    def __init__(self, row: str) -> None:
        parsed = self._parse_row(row)
        self._name = parsed["name"]
        self._index = int(parsed["index"])

    def _parse_row(self, row: str) -> str:
        patch = (
            Suppress("(") +
            Word(pp.printables).set_results_name("name") +
            Word(pp.nums).set_results_name("index") +
            Suppress(")")
        )
        return patch.parse_string(row).as_dict()

    def __repr__(self) -> str:
        return f"<Patch name=\"{self._name}\" index=\"{self._index}\" />"

    @property
    def name(self) -> str:
        """ Provides access to patch name. """
        return self._name

    @property
    def index(self) -> int:
        """ Provides access to patch index. """
        return self._index


class FluentSchemeHeader:
    """ Parse header row from a Scheme table file. """
    def __init__(self, row: str) -> None:
        self._names = self._parse_row(row)

    def _parse_row(self, row: str) -> list[str]:
        name = Word(pp.printables, exclude_chars="()")
        return np.ravel(name.search_string(row)).tolist()

    def __repr__(self) -> str:
        return f"<Header names=\"{' '.join(self._names)}\" />"

    @property
    def names(self) -> list[str]:
        """ Provides access to column names. """
        return self._names


class FluentSchemeTableRow:
    """ Parse data row from a Scheme table file. """
    def __init__(self, row: str) -> None:
        self._values = self._parse_row(row)

    def _parse_row(self, row: str) -> list[float]:
        values = OneOrMore(pp.common.sci_real)
        return np.ravel(values.search_string(row)).tolist()

    def __repr__(self) -> str:
        return f"<Data values=\"{' '.join(map(str, self._values))}\" />"

    @property
    def values(self) -> list[float]:
        """ Provides access to row data. """
        return self._values


class FluentDpmFile:
    """ Parse a sampled Scheme DPM file generated by Fluent. """
    def __init__(self, fname: str) -> None:
        self._patch = None
        self._header = None
        self._ncols = None
        self._data = []

        with open(fname) as fp:
            for k, row in enumerate(fp.readlines()):
                self._parse_row(k, row)

    def _get_patch(self, row: str) -> None:
        self._patch = FluentSchemePatch(row)

    def _get_header(self, row: str) -> None:
        self._header = FluentSchemeHeader(row).names

        # XXX: special case => when running with DEM collision
        # model, everything is fine, but when normal DPM is run
        # the `name` of injection (in outer parenthesis?) also
        # figures in header; try to understand better...
        if self._header[-1] == "name":
            self._header.pop()

        self._ncols = len(self._header)

    def _parse_row(self, k: int, row: str) -> None:
        match k:
            case 0:
                self._get_patch(row)
                return
            case 1:
                self._get_header(row)
                return

        data = FluentSchemeTableRow(row)

        if len(data.values) != self._ncols:
            raise ValueError(f"Malformed file at row {k+1}")

        self._data.append(data.values)

    def to_dataframe(self):
        """ Convert parsed data to a data frame for analysis. """
        return pd.DataFrame(self._data, columns=self._header)


def convert_xy_to_dict(
        fname: str | Path
        ) -> dict[int, dict[str, list[float]]]:
    """ Convert Ansys XY format to a dictionary structure.

    Parameters
    ----------
    fname: str | Path
        Path to ".xy" file to be converted by function.

    Returns
    -------
    dict[int, dict[str, list[float]]]
        A dictionary with each particle from pathlines transformed
        into a dictionary of coordinates X and Y.
    """
    with open(fname, "r") as fp:
        text = fp.read()

    start_index = 0
    block_index = 1
    all_blocks = {}

    while True:
        print(f"Extracting block {block_index}")

        # Markers for block start and end.
        start = f"((xy/key/label \"particle-{block_index}\")"
        end = ")"

        # Find index of next block start.
        start_index = text.find(start, start_index)

        # If none was found, then leave the loop.
        if start_index < 0:
            print(f"Block {block_index} not found, leaving now.")
            break

        # Start search after the end of start marker.
        start_index += len(start)

        # Find index of next block end.
        end_index = text.find(end, start_index)

        # If none was found, file is corrupted.
        if end_index < 0:
            raise Exception("No end index found, file is corrupted")

        # Parse the interval as a data frame for post-processing.
        block = StringIO(text[start_index:end_index])
        df = pd.read_csv(block, sep="\t", header=None)

        # Add data to main dictionary.
        all_blocks[block_index] = {}
        all_blocks[block_index]["X"] = df[0].to_list()
        all_blocks[block_index]["Y"] = df[1].to_list()

        # Set start indexes for next block lookup.
        start_index = end_index + len(end)
        block_index += 1

    return all_blocks


def load_dpm_table(fname: str) -> pd.DataFrame:
    """ Wrapper for simple loading of data table. """
    return FluentDpmFile(fname).to_dataframe()
