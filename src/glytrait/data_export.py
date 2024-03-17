"""Export data to files.

This module contains a function `export_all` that exports all data to files.
The data is provided as an Iterable of tuples (filename, data object).
This module will handle different data types by using the appropriate exporter.

Examples:
    >>> from glytrait.data_export import export_all
    >>> data_to_export = [
    ...     ("data.csv", pd.DataFrame(...)),
    ... ]
    >>> export_all(data_to_export, base_path="data")
"""

from collections.abc import Iterable, Callable
from pathlib import Path
from typing import Type, TypeVar, Any

import pandas as pd

from glytrait.formula import TraitFormula

__all__ = ["export_all"]


def export_all(to_export: Iterable[tuple[str, Any]], base_path: str) -> None:
    """Export all data to the directory given by `base_path`.

    Notes:
        `base_path` will not be created if it does not exist.

    Args:
        to_export: An iterable of tuples (filename, data object).
        base_path: The folder to export data into.
    """
    for filename, data in to_export:
        exporter = _fetch_exporter(data)
        filepath = str(Path(base_path) / filename)
        exporter(data, filepath)


ExportRegistryKey = tuple[Type, bool]
"""key: the object type, value: whether it is a list of objects"""

Exporter = Callable[[Any, str], None]
"""A function that exports an object to a file.
The first argument is the object to export, the second is the filepath.
"""

export_registry: dict[ExportRegistryKey, Exporter] = {}


def register_exporter(data_type: Type, is_list: bool = False):  # type: ignore
    """Decorator for regitering an exporter for a certain type."""
    T = TypeVar("T", bound=Exporter)

    def decorator(exporter_class: T) -> T:
        export_registry[(data_type, is_list)] = exporter_class
        return exporter_class

    return decorator


def _fetch_exporter(single_data: Any) -> Exporter:
    """Factory function for creating an exporter for a data object."""
    is_list = isinstance(single_data, list)
    if is_list:
        try:
            data_type = type(single_data[0])
        except IndexError:
            return dummy_exporter
    else:
        data_type = type(single_data)
    exporter = export_registry.get((data_type, is_list))
    if exporter is not None:
        return exporter
    else:
        raise ValueError(f"Unsupported data type for export: {data_type}.")


def dummy_exporter(data: Any, filepath: str) -> None:
    """Dummy exporter that does nothing."""
    pass


@register_exporter(pd.DataFrame)
def export_dataframe(df: pd.DataFrame, filepath: str) -> None:
    """Export a pandas DataFrame to a CSV file."""
    df.to_csv(filepath, index=True)
