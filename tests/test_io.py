import os
from unittest.mock import MagicMock, patch, call

import numpy as np
import pandas as pd
import pytest
from openpyxl import load_workbook

from glytrait import glycan as glyc
from glytrait import io
from glytrait import trait
from glytrait.exception import *


def test_read_input(mocker):
    glycoct_parser_mock = mocker.patch(
        "glytrait.glycan.NGlycan.from_glycoct",
        return_value=MagicMock(spec=glyc.NGlycan),
    )
    test_data = {
        "Glycan ID": ["H3N4S1", "H5N4S1"],
        "Structure": ["some_glycoct", "another_glycoct"],
        "Sample1": [2.3, 1.2],
        "Sample2": [1.2, 2.3],
        "Sample3": [3.4, 3.4],
    }
    test_df = pd.DataFrame(test_data)
    expected_abundance_table = (
        test_df.drop(columns=["Structure"]).set_index("Glycan ID").T
    )

    with patch("pandas.read_csv", return_value=test_df.copy()):
        glycans, abundance_table = io.read_input("fake_path")

    glycoct_parser_mock.assert_has_calls(
        [call(structure) for structure in test_df["Structure"]]
    )
    assert len(glycans) == len(test_df)
    assert all(isinstance(g, glyc.NGlycan) for g in glycans)
    pd.testing.assert_frame_equal(abundance_table, expected_abundance_table)


def test_read_input_with_structure_parse_error():
    test_data = {
        "Glycan ID": ["H3N4S1", "H5N4S1"],
        "Structure": ["some_glycoct", "another_glycoct"],
        "Sample1": [2.3, 1.2],
        "Sample2": [1.2, 2.3],
        "Sample3": [3.4, 3.4],
    }
    test_df = pd.DataFrame(test_data)

    with patch("pandas.read_csv", return_value=test_df.copy()):
        with pytest.raises(StructureParseError) as excinfo:
            io.read_input("fake_path")
        assert "H3N4S1" in str(excinfo.value)


@pytest.mark.parametrize(
    "invalid_df, error_message", [
        (pd.DataFrame({"Invalid ID": ["glycan1", "glycan2"], "Structure": ["structure1", "structure2"]}),
         "The first column of the input file should be Glycan ID."),
        (pd.DataFrame({"Glycan ID": ["glycan1", "glycan2"], "Invalid Structure": ["structure1", "structure2"]}),
         "The second column of the input file should be Structure."),
        (pd.DataFrame({"Glycan ID": ["glycan1", "glycan1"], "Structure": ["structure1", "structure2"]}),
         "There are duplicated glycan IDs in the input file."),
        (pd.DataFrame({"Glycan ID": ["glycan1", "glycan2"], "Structure": ["structure1", "structure1"]}),
         "There are duplicated structures in the input file."),
        (pd.DataFrame({"Glycan ID": ["glycan1", "glycan2"], "Structure": ["structure1", "structure2"], "Sample1": ["invalid", "data"]}),
         "The abundance columns in the input file should be numeric."),
        (pd.DataFrame({"Glycan ID": ["glycan1", np.nan], "Structure": ["structure1", "structure2"]}),
         "There are missing values in the input file.")
    ]
)
def test_read_input_with_invalid_data(invalid_df, error_message, clean_dir):
    file = clean_dir / "invalid_data.csv"
    invalid_df.to_csv(str(file), index=False)
    with pytest.raises(InputError, match=error_message):
        io.read_input(str(file))
