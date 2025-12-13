from typing import Any
from dataclasses import dataclass, field

__all__ = [
    "LookupTable"
]


@dataclass(frozen=True)
class LookupTable:
    """
    Implements a lookup table with a column header, a row header and a data
    table.

    Two lookup mechanisms are provided:
    1.  Given a column header value and a row header value, returns the value
        in the data cell whose column index and row index correspond with the
        index of the column header value and the index of the row header value.

    2.  Given a column header value and a data value, returns the value in the
        row header whose index is equal to the row index of the selected
        cell in the data table. The column index of the selected cell in the
        data table is equal to the index of the given column header value. The
        row index of the selected cell is determined by taking the index of the
        data value in the selected column that is equal or just greater than the
        given data value. Note that data values in the data table columns must
        be in ascending order.
    """
    data: list[list[float]] = field(default_factory=list)
    col_header: list[Any] | dict[Any, str] = field(default_factory=list)
    row_header: list[Any] | dict[Any, str] = field(default_factory=list)
    description: str = ""
    cols_description: str = ""
    rows_description: str = ""

    @classmethod
    def create(
        cls,
        row_header: list[Any] | dict[Any, str],
        col_header: list[Any] | dict[Any, str],
        data: list[list[float]],
        description: str = "",
        rows_description: str = "",
        cols_description: str = ""
    ) -> 'LookupTable':
        """
        Creates a lookup table.

        Parameters
        ----------
        row_header: list[Any] | dict[Any, str]
            List with the row header values or a dict of which the keys are
            the row header values and the dict values are strings that explain
            the meaning of the row header values.
        col_header: list[Any] | dict[Any, str]
            List the column header values or a dict of which the keys are
            the column header values and the dict values are strings that
            explain the meaning of the column header values.
        data: list[list[float]]
            A table of (floating) numbers implemented as a lists in a list.
        description: str, optional
            Describes the contents of the lookup table.
        rows_description: str, optional
            Describes the meaning of the row header values.
        cols_description: str, optional
            Describes the meaning of the column header values.

        Returns
        -------
        LookupTable
        """
        if len(data) != len(row_header):
            raise ValueError(
                f"The number of rows in `data` ({len(data)}) is not equal "
                f"to the length of `row_header` ({len(row_header)})."
            )
        for i, row in enumerate(data):
            if len(row) < len(col_header):
                d = len(col_header) - len(row)
                row.extend([float("nan")] * d)
            if len(row) > len(col_header):
                raise ValueError(
                    f"The number of values ({len(row)}) in row {i} is greater "
                    f"than the number of columns ({len(col_header)})."
                )
        return cls(
            data,
            col_header,
            row_header,
            description,
            cols_description,
            rows_description
        )

    def data_value(self, rowheader_val: Any, colheader_val: Any) -> float:
        """
        Returns the data value that belongs to the given row header value and
        the given column header value.

        Parameters
        ----------
        rowheader_val: Any
            Value present in the row header of the lookup table.
        colheader_val: Any
            Value present in the column header of the lookup table.

        Returns
        -------
        float
        """
        if isinstance(self.row_header, dict):
            row_header = list(self.row_header.keys())
        else:
            row_header = self.row_header
        if isinstance(self.col_header, dict):
            col_header = list(self.col_header.keys())
        else:
            col_header = self.col_header
        i = row_header.index(rowheader_val)
        j = col_header.index(colheader_val)
        value = self.data[i][j]
        return value

    def rowheader_value(self, colheader_val: Any, data_val: float) -> Any:
        """
        Returns the row header value that belongs to the given column header
        value and the given data value.

        To determine the row header value, first the column of the data table is
        selected that has the same index as the index of the given column header
        value. Then the index is determined of the data value in this column
        which is equal to or just greater than the given data value. This index
        finally gives us the row header value. Data values in the data table
        columns must be in ascending order.

        Parameters
        ----------
        colheader_val: Any
            Value present in the column header of the lookup table.
        data_val: float
            Value for which we want the corresponding row header value.

        Returns
        -------
        Any
        """
        if isinstance(self.col_header, dict):
            col_header = list(self.col_header.keys())
        else:
            col_header = self.col_header
        j = col_header.index(colheader_val)
        data_vals = [self.data[i][j] for i in range(len(self.data))]
        for i in range(len(data_vals)):
            if i < len(data_vals) - 1:
                if data_vals[i] < data_val <= data_vals[i + 1]:
                    if isinstance(self.row_header, dict):
                        row_header = list(self.row_header.keys())
                        return row_header[i + 1]
                    else:
                        return self.row_header[i + 1]
        else:
            raise ValueError(
                "The given data value falls outside the data table"
            )

    def colheader_value(self, rowheader_val: Any, data_val: float) -> Any:
        """
        Returns the column header value that belongs to the given row header
        value and the given data value.

        To determine the column header value, first the row of the data table is
        selected that has the same index as the index of the given row header
        value. Then the index is determined of the data value in this row which
        is equal to or just greater than the given data value. This index
        finally gives us the column header value. Data values in the data table
        rows must be in ascending order.

        Parameters
        ----------
        rowheader_val: Any
            Value present in the row header of the lookup table.
        data_val: float
            Value for which we want the corresponding column header value.

        Returns
        -------
        Any
        """
        if isinstance(self.row_header, dict):
            row_header = list(self.row_header.keys())
        else:
            row_header = self.row_header
        i = row_header.index(rowheader_val)
        data_vals = self.data[i]
        for j in range(len(data_vals)):
            if j < len(data_vals) - 1:
                if data_vals[j] < data_val <= data_vals[j + 1]:
                    if isinstance(self.col_header, dict):
                        col_header = list(self.col_header.keys())
                        return col_header[j + 1]
                    else:
                        return self.col_header[j + 1]
        else:
            raise ValueError(
                "The given data value falls outside the data table"
            )
