import pandas as pd
import pytest

from glytrait import io
from glytrait.exception import *
from glytrait.io import load_default_structures

from .glycoct import *

pytestmark = pytest.mark.skip("Not implemented yet.")


def test_read_group_file(mocker):
    group_df = pd.DataFrame(
        {
            "Group": ["Group1", "Group2", "Group1", "Group2"],
        },
        index=["Sample1", "Sample2", "Sample3", "Sample4"],
    )
    expected = pd.Series(group_df["Group"], index=group_df.index)
    mocker.patch("pandas.read_csv", return_value=group_df)
    pd.testing.assert_series_equal(io.read_group_file("fake_path"), expected)


def test_read_group_file_from_df():
    group_df = pd.DataFrame(
        {
            "Group": ["Group1", "Group2", "Group1", "Group2"],
        },
        index=["Sample1", "Sample2", "Sample3", "Sample4"],
    )
    expected = pd.Series(group_df["Group"], index=group_df.index)
    pd.testing.assert_series_equal(io.read_group_file(group_df), expected)


def test_read_group_file_wrong_format(mocker):
    group_df = pd.DataFrame(
        {
            "Group": ["Group1", "Group2", "Group1", "Group2"],
            "Wrong": ["Wrong", "Wrong", "Wrong", "Wrong"],
        },
        index=["Sample1", "Sample2", "Sample3", "Sample4"],
    )
    mocker.patch("pandas.read_csv", return_value=group_df)
    with pytest.raises(InputError) as excinfo:
        io.read_group_file("fake_path")
    assert "The group file should only have two columns." in str(excinfo.value)


def test_read_structure_file_from_csv(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    result = io.read_structure_file("fake_path", ["glycan2", "glycan1"])
    assert len(result) == 2


def test_read_structure_file_from_df():
    structure_df = pd.DataFrame(
        {
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    result = io.read_structure_file(structure_df, ["glycan2", "glycan1"])
    assert len(result) == 2


def test_read_structure_file_from_dir(mocker, clean_dir):
    glycan_1 = clean_dir / "glycan1.glycoct_condensed"
    glycan_2 = clean_dir / "glycan2.glycoct_condensed"
    glycan_3 = clean_dir / "glycan3.glycoct_condensed"

    glycan_1.write_text(test_glycoct_1)
    glycan_2.write_text(test_glycoct_2)
    glycan_3.write_text(test_glycoct_3)

    result = io.read_structure_file(clean_dir, ["glycan2", "glycan1"])
    assert len(result) == 2


def test_read_structure_file_wrong_format(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
            "Wrong": ["Wrong", "Wrong", "Wrong"],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    with pytest.raises(InputError) as excinfo:
        io.read_structure_file("fake_path", ["glycan2", "glycan1"])
    assert "The structure file should only have two columns." in str(excinfo.value)


def test_read_structure_file_with_missing_structure(mocker):
    structure_df = pd.DataFrame(
        {
            "Structure": [test_glycoct_1, test_glycoct_2, test_glycoct_3],
        },
        index=["glycan1", "glycan2", "glycan3"],
    )
    mocker.patch("pandas.read_csv", return_value=structure_df)
    with pytest.raises(InputError) as excinfo:
        io.read_structure_file("fake_path", ["glycan2", "glycan1", "glycan4"])
    assert "glycan4 is not found in the structure file." in str(excinfo.value)


@pytest.mark.parametrize("db", ["serum", "IgG"])
def test_load_default(db):
    result = load_default_structures(db, ["H3N3", "H5N4S2"])
    assert len(result) == 2
