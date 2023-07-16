import pandas as pd
import pytest

from glytrait import io
from glytrait.exception import *
from glytrait.io import load_default_structures

from .glycoct import *


def test_read_input_with_structure_col(mocker):
    df = pd.DataFrame(
        {
            "Composition": ["comp1", "comp2", "comp3"],
            "Structure": ["structure1", "structure2", "structure3"],
            "Sample1": [1., 2., 3.],
            "Sample2": [4., 5., 6.],
            "Sample3": [7., 8., 9.],
        }
    )
    expected = (
        ["comp1", "comp2", "comp3"],
        ["structure1", "structure2", "structure3"],
        pd.DataFrame(
            {
                "comp1": [1., 4., 7.],
                "comp2": [2., 5., 8.],
                "comp3": [3., 6., 9.],
            },
            index=["Sample1", "Sample2", "Sample3"],
        ),
    )
    mocker.patch("pandas.read_csv", return_value=df)
    result = io.read_input("fake_path")
    assert result[0] == expected[0]
    assert result[1] == expected[1]
    pd.testing.assert_frame_equal(result[2], expected[2], check_names=False)


def test_read_input_without_structure_col(mocker):
    df = pd.DataFrame(
        {
            "Composition": ["comp1", "comp2", "comp3"],
            "Sample1": [1., 2., 3.],
            "Sample2": [4., 5., 6.],
            "Sample3": [7., 8., 9.],
        }
    )
    expected = (
        ["comp1", "comp2", "comp3"],
        None,
        pd.DataFrame(
            {
                "comp1": [1., 4., 7.],
                "comp2": [2., 5., 8.],
                "comp3": [3., 6., 9.],
            },
            index=["Sample1", "Sample2", "Sample3"],
        ),
    )
    mocker.patch("pandas.read_csv", return_value=df)
    result = io.read_input("fake_path")
    assert result[0] == expected[0]
    assert result[1] == expected[1]
    pd.testing.assert_frame_equal(result[2], expected[2], check_names=False)


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
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    result = io.read_structure("fake_path", ["glycan2", "glycan1"])
    assert len(result) == 2


def test_read_structure_wrong_format(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
            "Wrong": ["Wrong", "Wrong", "Wrong"],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    with pytest.raises(InputError) as excinfo:
        io.read_structure("fake_path", ["glycan2", "glycan1"])
    assert "The structure file should only have two columns." in str(excinfo.value)


def test_read_structure_with_missing_structure(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    with pytest.raises(InputError) as excinfo:
        io.read_structure("fake_path", ["glycan2", "glycan1", "glycan4"])
    assert "glycan4 is not found in the structure file." in str(excinfo.value)


@pytest.mark.parametrize("db", ["serum", "IgG"])
def test_load_default(db):
    result = load_default_structures(db, ["H3N3", "H5N4S2"])
    assert len(result) == 2
