from typing import NoReturn, Optional

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Side
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.worksheet.worksheet import Worksheet

from glytrait.config import Config
from glytrait.exception import InputError, StructureParseError
from glytrait.glycan import NGlycan
from glytrait.trait import TraitFormula


def read_input(file: str) -> tuple[list[NGlycan], pd.DataFrame]:
    """Read the input file.

    The first column should be the glycan IDs. The composotion will do. The second column should
    be the glycan structure strings. The rest columns should be the abundance of the glycans in
    different samples.

    An example of the input file:
    ```
    Glycan ID, Structure, Sample1, Sample2, Sample3
    H3N4S1, some_glycoct, 2.3, 1.2, 3.4
    H5N4S1, another_glycoct, 1.2, 2.3, 3.4
    ```

    Args:
        file (str): The input csv or xlsx file.

    Returns:
        glycans (list[NGlycan]): The glycans.
        abundance_table (pd.DataFrame): The abundance table, with samples as index and glycan IDs
            as columns.

    Raises:
        InputError: If the input file format is wrong.
    """
    df = pd.read_csv(file)
    df = df.dropna(how="all")
    _check_input(df)
    df = df.set_index("Glycan ID")
    glycans: list[NGlycan] = []
    for glycan_id, structure in zip(df.index, df["Structure"]):
        try:
            glycan = NGlycan.from_glycoct(structure)
        except StructureParseError as e:
            raise StructureParseError(glycan_id + ": " + str(e))
        else:
            glycans.append(glycan)
    abundance_table = df.drop(columns=["Structure"]).T
    return glycans, abundance_table


def read_group(file: str) -> pd.Series:
    """Read the group file.

    The first column should be the sample IDs and the second column should be the group names.

    Args:
        file (str): The group file.

    Returns:
        groups (pd.Series): The group names, with sample IDs as index.
    """
    df = pd.read_csv(file, index_col=0)
    if df.shape[1] != 1:
        raise InputError("The group file should only have two columns.")
    groups = df.squeeze()
    return groups


def _check_input(df: pd.DataFrame) -> NoReturn:
    """Check the input dataframe."""
    if df.columns[0] != "Glycan ID":
        raise InputError("The first column of the input file should be Glycan ID.")
    if df.columns[1] != "Structure":
        raise InputError("The second column of the input file should be Structure.")
    if df["Glycan ID"].duplicated().sum() > 0:
        raise InputError("There are duplicated glycan IDs in the input file.")
    if df["Structure"].duplicated().sum() > 0:
        raise InputError("There are duplicated structures in the input file.")
    if np.any(df.iloc[:, 2:].dtypes != "float64"):
        raise InputError("The abundance columns in the input file should be numeric.")
    if df.isna().sum().sum() > 0:
        raise InputError("There are missing values in the input file.")


def write_output(
    config: Config,
    derived_traits: pd.DataFrame,
    direct_traits: pd.DataFrame,
    meta_prop_df: pd.DataFrame,
    formulas: list[TraitFormula],
    hypothesis_test: Optional[pd.DataFrame] = None,
) -> None:
    """Write the output file.

    Args:
        config (Config): The configuration.
        derived_traits (pd.DataFrame): The trait values, with samples as index and trait names as
            columns.
        direct_traits (pd.DataFrame): The normalized abundance table, with samples as index and
            glycan IDs as columns.
        meta_prop_df (pd.DataFrame): The table of meta properties generated by
            `build_meta_property_table`.
        formulas (list[TraitFormula]): The trait formulas.
        hypothesis_test (Optional[pd.DataFrame], optional): The hypothesis test results. Defaults
            to None.
    """
    wb = Workbook()

    # The abundance table
    ws1: Worksheet = wb.active
    ws1.title = "Trait values"
    combined_df = pd.concat([direct_traits, derived_traits], axis=1)

    for row in dataframe_to_rows(combined_df, index=True, header=True):
        ws1.append(row)
    ws1.delete_rows(2)
    for cell in ws1[1][1:]:
        cell.font = cell.font.copy(bold=True)
        cell.border = cell.border.copy(bottom=Side(border_style="double"))
    ws1.insert_rows(1)
    ws1.merge_cells(
        start_row=1,
        start_column=2,
        end_row=1,
        end_column=len(direct_traits.columns) + 1,
    )
    ws1.cell(1, 2).value = "Direct traits"
    ws1.merge_cells(
        start_row=1,
        start_column=len(direct_traits.columns) + 2,
        end_row=1,
        end_column=len(combined_df.columns) + 1,
    )
    ws1.cell(1, len(direct_traits.columns) + 2).value = "Derived traits"
    for cell in ws1[1][1:]:
        cell.font = cell.font.copy(bold=True)
        cell.border = cell.border.copy(bottom=Side(border_style="thin"))

    # The trait definitions
    ws2 = wb.create_sheet("Trait definitions")
    ws2.append(["Trait Name", "Description"])
    for cell in ws2[1]:
        cell.font = cell.font.copy(bold=True)
    for formula in formulas:
        ws2.append([formula.name, formula.description])

    # The meta properties
    ws3 = wb.create_sheet("Meta properties")
    for row in dataframe_to_rows(meta_prop_df, index=True, header=True):
        ws3.append(row)
    ws3.delete_rows(2)
    for cell in ws3[1]:
        cell.font = cell.font.copy(bold=True)

    # The hypothesis test results
    if hypothesis_test is not None:
        ws4 = wb.create_sheet("Hypothesis test results")
        for row in dataframe_to_rows(hypothesis_test, index=True, header=True):
            ws4.append(row)
        ws4.delete_rows(2)
        for cell in ws4[1]:
            cell.font = cell.font.copy(bold=True)

    wb.save(config.get("output_file"))
