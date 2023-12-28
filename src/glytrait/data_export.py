from collections.abc import Iterable, Callable
from pathlib import Path
from typing import Type, TypeVar, Any

import pandas as pd

from glytrait.formula import TraitFormula


def export_all(to_export: Iterable[tuple[str, Any]], base_path: str) -> None:
    """Export all data to the directory given by `base_path`.

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


@register_exporter(TraitFormula, is_list=True)
def export_formulas(formulas: list[TraitFormula], filepath: str) -> None:
    """Export a list of TraitFormulas to a CSV file."""
    with open(filepath, "w", encoding="utf8") as f:
        for formula in formulas:
            row = f"{formula.name}: {formula.description}\n"
            f.write(row)
