from unittest.mock import MagicMock, patch, call

import numpy as np
import pandas as pd
import pytest

from glytrait import glycan as glyc
from glytrait import io
from glytrait.exception import *


def test_read_input_with_structure(mocker):
    glycoct_parser_mock = mocker.patch(
        "glytrait.glycan.NGlycan.from_glycoct",
        return_value=MagicMock(spec=glyc.NGlycan),
    )
    test_data = {
        "Composition": ["H3N4S1", "H5N4S1"],
        "Structure": ["some_glycoct", "another_glycoct"],
        "Sample1": [2.3, 1.2],
        "Sample2": [1.2, 2.3],
        "Sample3": [3.4, 3.4],
    }
    test_df = pd.DataFrame(test_data)
    expected_abundance_table = (
        test_df.drop(columns=["Structure"]).set_index("Composition").T
    )

    with patch("pandas.read_csv", return_value=test_df.copy()):
        glycans, abundance_table = io.read_input("fake_path")

    glycoct_parser_mock.assert_has_calls(
        [call(structure) for structure in test_df["Structure"]]
    )
    assert len(glycans) == len(test_df)
    assert all(isinstance(g, glyc.NGlycan) for g in glycans)
    pd.testing.assert_frame_equal(abundance_table, expected_abundance_table)


def test_read_input_without_structure():
    test_data = {
        "Composition": ["H3N4S1", "H5N4S1"],
        "Sample1": [2.3, 1.2],
        "Sample2": [1.2, 2.3],
        "Sample3": [3.4, 3.4],
    }
    test_df = pd.DataFrame(test_data)
    expected_abundance_table = test_df.set_index("Composition").T

    with patch("pandas.read_csv", return_value=test_df.copy()):
        glycans, abundance_table = io.read_input("fake_path")

    assert glycans is None
    pd.testing.assert_frame_equal(abundance_table, expected_abundance_table)


def test_read_input_with_structure_parse_error():
    test_data = {
        "Composition": ["H3N4S1", "H5N4S1"],
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
    "invalid_df, error_message",
    [
        (
            pd.DataFrame(
                {
                    "Invalid ID": ["glycan1", "glycan2"],
                    "Structure": ["structure1", "structure2"],
                }
            ),
            "The first column of the input file should be Composition.",
        ),
        (
            pd.DataFrame(
                {
                    "Composition": ["glycan1", "glycan2"],
                    "Invalid Structure": ["structure1", "structure2"],
                }
            ),
            "The abundance columns in the input file should be numeric.",
        ),
        (
            pd.DataFrame(
                {
                    "Composition": ["glycan1", "glycan1"],
                    "Structure": ["structure1", "structure2"],
                }
            ),
            "There are duplicated Compositions in the input file.",
        ),
        (
            pd.DataFrame(
                {
                    "Composition": ["glycan1", "glycan2"],
                    "Structure": ["structure1", "structure1"],
                }
            ),
            "There are duplicated Structures in the input file.",
        ),
        (
            pd.DataFrame(
                {
                    "Composition": ["glycan1", "glycan2"],
                    "Structure": ["structure1", "structure2"],
                    "Sample1": ["invalid", "data"],
                }
            ),
            "The abundance columns in the input file should be numeric.",
        ),
        (
            pd.DataFrame(
                {
                    "Composition": ["glycan1", np.nan],
                    "Structure": ["structure1", "structure2"],
                }
            ),
            "There are missing values in the input file.",
        ),
    ],
)
def test_read_input_with_invalid_data(invalid_df, error_message, clean_dir):
    file = clean_dir / "invalid_data.csv"
    invalid_df.to_csv(str(file), index=False)
    with pytest.raises(InputError, match=error_message):
        io.read_input(str(file))


def test_read_group(mocker):
    group_df = pd.DataFrame(
        {
            "Group": ["Group1", "Group2", "Group1", "Group2"],
        },
        index=["Sample1", "Sample2", "Sample3", "Sample4"],
    )
    expected = pd.Series(group_df["Group"], index=group_df.index)
    mocker.patch("pandas.read_csv", return_value=group_df)
    pd.testing.assert_series_equal(io.read_group("fake_path"), expected)


def test_read_group_wrong_format(mocker):
    group_df = pd.DataFrame(
        {
            "Group": ["Group1", "Group2", "Group1", "Group2"],
            "Wrong": ["Wrong", "Wrong", "Wrong", "Wrong"],
        },
        index=["Sample1", "Sample2", "Sample3", "Sample4"],
    )
    mocker.patch("pandas.read_csv", return_value=group_df)
    with pytest.raises(InputError) as excinfo:
        io.read_group("fake_path")
    assert "The group file should only have two columns." in str(excinfo.value)


def test_read_structure(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": ["structure1", "structure2", "structure3"],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    expected = ["structure2", "structure1"]

    mocker.patch("pandas.read_csv", return_value=structure_df)
    mocker.patch("glytrait.io.NGlycan.from_glycoct", side_effect=lambda x: x)
    result = io.read_structure("fake_path", ["glycan2", "glycan1"])
    assert result == expected


def test_read_structure_wrong_format(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": ["structure1", "structure2", "structure3"],
            "Wrong": ["Wrong", "Wrong", "Wrong"],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    mocker.patch("glytrait.io.NGlycan.from_glycoct", side_effect=lambda x: x)
    with pytest.raises(InputError) as excinfo:
        io.read_structure("fake_path", ["glycan2", "glycan1"])
    assert "The structure file should only have two columns." in str(excinfo.value)


def test_read_structure_with_missing_structure(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": ["structure1", "structure2", "structure3"],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    mocker.patch("glytrait.io.NGlycan.from_glycoct", side_effect=lambda x: x)
    with pytest.raises(InputError) as excinfo:
        io.read_structure("fake_path", ["glycan2", "glycan1", "glycan4"])
    assert "glycan4 is not found in the structure file." in str(excinfo.value)
