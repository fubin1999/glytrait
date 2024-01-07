from pathlib import Path

import pandas as pd
import pytest

from glytrait.exception import (
    FileTypeError,
    FileFormatError,
    StructureParseError,
    NotEnoughGroupsError,
)
from glytrait.load_data import (
    load_input_data,
    load_abundance_table,
    load_structures,
    load_compositions,
    load_groups,
    GlyTraitInputData,
)


@pytest.fixture
def abundance_table() -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "G1": [1, 2, 3],
            "G2": [4, 5, 6],
            "G3": [7, 8, 9],
        },
        index=pd.Index(["sample1", "sample2", "sample3"], name="Sample"),
    )
    return df


@pytest.fixture
def abundance_table_file(abundance_table, clean_dir) -> str:
    filepath = str(clean_dir / "abundance_table.csv")
    abundance_table.to_csv(filepath, index=True)
    return filepath


class TestLoadAbundanceTable:
    def test_basic(self, abundance_table_file, abundance_table):
        result = load_abundance_table(abundance_table_file)
        pd.testing.assert_frame_equal(result, abundance_table)

    def test_not_csv(self, clean_dir):
        filepath = clean_dir / "abundance_table.txt"
        filepath.touch()
        with pytest.raises(FileTypeError):
            load_abundance_table(filepath)

    def test_not_exist(self, clean_dir):
        filepath = clean_dir / "abundance_table.csv"
        with pytest.raises(FileNotFoundError):
            load_abundance_table(filepath)

    def test_not_numeric(self, abundance_table, clean_dir):
        abundance_table["sample1"] = ["a", "b", "c"]
        filepath = str(clean_dir / "abundance_table.csv")
        abundance_table.to_csv(filepath, index=True)
        with pytest.raises(FileFormatError):
            load_abundance_table(filepath)

    def test_duplicated_glycans(self, abundance_table, clean_dir):
        abundance_table.index = pd.Index(["G1", "G1", "G3"], name="GlycanID")
        filepath = str(clean_dir / "abundance_table.csv")
        abundance_table.to_csv(filepath, index=True)
        with pytest.raises(FileFormatError):
            load_abundance_table(filepath)

    def test_duplicated_sample_names(self, abundance_table, clean_dir):
        abundance_table.columns = ["sample1", "sample1", "sample3"]
        filepath = str(clean_dir / "abundance_table.csv")
        abundance_table.to_csv(filepath, index=True)
        with pytest.raises(FileFormatError):
            load_abundance_table(filepath)

    def test_no_glycan_id_column(self, abundance_table, clean_dir):
        abundance_table.index.name = "GlycanID1"
        filepath = str(clean_dir / "abundance_table.csv")
        abundance_table.to_csv(filepath, index=True)
        with pytest.raises(FileFormatError):
            load_abundance_table(filepath)


@pytest.fixture
def structure_file(clean_dir) -> str:
    filepath = str(clean_dir / "structure.csv")
    df = pd.DataFrame(
        {
            "GlycanID": ["G1", "G2", "G3"],
            "Structure": ["glycoct1", "glycoct2", "glycoct3"],
        }
    )
    df.to_csv(filepath, index=False)
    return filepath


@pytest.fixture
def patch_structure_from_string(monkeypatch):
    def mock_return(name, string):
        if string == "invalid":
            raise StructureParseError
        return string

    monkeypatch.setattr("glytrait.glycan.Structure.from_string", mock_return)


@pytest.fixture
def structure_folder(clean_dir) -> str:
    dirpath = clean_dir / "structures"
    dirpath.mkdir()
    for i in range(1, 4):
        filepath = dirpath / f"G{i}.glycoct_condensed"
        filepath.write_text(f"glycoct{i}")
    return str(dirpath)


@pytest.mark.usefixtures("patch_structure_from_string")
class TestLoadStructures:
    def test_basic_from_csv(self, structure_file):
        result = load_structures(structure_file)
        assert result == {
            "G1": "glycoct1",
            "G2": "glycoct2",
            "G3": "glycoct3",
        }

    def test_not_csv(self, clean_dir):
        filepath = clean_dir / "structure.txt"
        filepath.touch()
        with pytest.raises(FileTypeError):
            load_structures(filepath)

    def test_not_exist(self, clean_dir):
        filepath = clean_dir / "structure.csv"
        with pytest.raises(FileNotFoundError):
            load_structures(filepath)

    @pytest.mark.parametrize("column", ["GlycanID", "Structure"])
    def test_missing_columns(self, structure_file, clean_dir, column):
        df = pd.read_csv(structure_file)
        df = df.drop(column, axis=1)
        filepath = str(clean_dir / "structure.csv")
        df.to_csv(filepath, index=False)
        with pytest.raises(FileFormatError):
            load_structures(filepath)

    def test_duplicated_glycan_ids(self, structure_file, clean_dir):
        df = pd.read_csv(structure_file)
        df.loc[2, "GlycanID"] = "G1"
        filepath = str(clean_dir / "structure.csv")
        df.to_csv(filepath, index=False)
        with pytest.raises(FileFormatError):
            load_structures(filepath)

    def test_duplicated_structure_strings(self, structure_file, clean_dir):
        df = pd.read_csv(structure_file)
        df.loc[2, "Structure"] = "glycoct1"
        filepath = str(clean_dir / "structure.csv")
        df.to_csv(filepath, index=False)
        with pytest.raises(FileFormatError):
            load_structures(filepath)

    def test_from_dir(self, structure_folder):
        result = load_structures(structure_folder)
        assert result == {
            "G1": "glycoct1",
            "G2": "glycoct2",
            "G3": "glycoct3",
        }

    def test_from_dir_not_exist(self, clean_dir):
        filepath = clean_dir / "structures"
        with pytest.raises(FileNotFoundError):
            load_structures(filepath)

    def test_from_dir_duplicate_structures(self, structure_folder):
        filepath = structure_folder + "/G1.glycoct_condensed"
        Path(filepath).write_text("glycoct2")
        with pytest.raises(FileFormatError):
            load_structures(structure_folder)


@pytest.fixture
def composition_file(clean_dir) -> str:
    filepath = str(clean_dir / "composition.csv")
    df = pd.DataFrame(
        {
            "GlycanID": ["G1", "G2", "G3"],
            "Composition": ["comp1", "comp2", "comp3"],
        }
    )
    df.to_csv(filepath, index=False)
    return filepath


@pytest.fixture
def patch_composition_from_string(monkeypatch):
    def mock_return(name, string):
        if string == "invalid":
            raise StructureParseError
        return string

    monkeypatch.setattr("glytrait.glycan.Composition.from_string", mock_return)


@pytest.mark.usefixtures("patch_composition_from_string")
class TestLoadCompositions:
    def test_basic(self, composition_file):
        result = load_compositions(composition_file)
        assert result == {
            "G1": "comp1",
            "G2": "comp2",
            "G3": "comp3",
        }

    def test_not_csv(self, clean_dir):
        filepath = clean_dir / "composition.txt"
        filepath.touch()
        with pytest.raises(FileTypeError):
            load_compositions(filepath)

    def test_not_exist(self, clean_dir):
        filepath = clean_dir / "composition.csv"
        with pytest.raises(FileNotFoundError):
            load_compositions(filepath)

    @pytest.mark.parametrize("column", ["GlycanID", "Composition"])
    def test_missing_columns(self, composition_file, clean_dir, column):
        df = pd.read_csv(composition_file)
        df = df.drop(column, axis=1)
        filepath = str(clean_dir / "composition.csv")
        df.to_csv(filepath, index=False)
        with pytest.raises(FileFormatError):
            load_compositions(filepath)

    def test_duplicated_glycan_ids(self, composition_file, clean_dir):
        df = pd.read_csv(composition_file)
        df.loc[2, "GlycanID"] = "G1"
        filepath = str(clean_dir / "composition.csv")
        df.to_csv(filepath, index=False)
        with pytest.raises(FileFormatError):
            load_compositions(filepath)

    def test_duplicated_compositions(self, composition_file, clean_dir):
        df = pd.read_csv(composition_file)
        df.loc[2, "Composition"] = "comp1"
        filepath = str(clean_dir / "composition.csv")
        df.to_csv(filepath, index=False)
        with pytest.raises(FileFormatError):
            load_compositions(filepath)


@pytest.fixture
def group_series() -> pd.DataFrame:
    s = pd.Series(
        ["group1", "group2", "group3"],
        name="Group",
        index=pd.Index(["sample1", "sample2", "sample3"], name="Sample"),
    )
    return s


@pytest.fixture
def group_series_file(group_series, clean_dir) -> str:
    filepath = str(clean_dir / "group_table.csv")
    group_series.to_csv(filepath, index=True)
    return filepath


class TestLoadGroupTable:
    def test_basic(self, group_series_file, group_series):
        result = load_groups(group_series_file)
        pd.testing.assert_series_equal(result, group_series)

    def test_not_csv(self, clean_dir):
        filepath = clean_dir / "group_table.txt"
        filepath.touch()
        with pytest.raises(FileTypeError):
            load_groups(filepath)

    def test_not_exist(self, clean_dir):
        filepath = clean_dir / "group_table.csv"
        with pytest.raises(FileNotFoundError):
            load_groups(filepath)

    def test_duplicated_samples(self, group_series, clean_dir):
        group_series.index = pd.Index(["sample1", "sample1", "sample3"], name="Sample")
        filepath = str(clean_dir / "group_table.csv")
        group_series.to_csv(filepath, index=True)
        with pytest.raises(FileFormatError):
            load_groups(filepath)

    @pytest.mark.parametrize("column", ["Group", "Sample"])
    def test_missing_columns(self, group_series, clean_dir, column):
        df = group_series.reset_index()
        df = df.drop(column, axis=1)
        filepath = str(clean_dir / "group_table.csv")
        df.to_csv(filepath, index=True)
        with pytest.raises(FileFormatError):
            load_groups(filepath)

    def test_group_number(self, group_series, clean_dir):
        group_series = pd.Series(
            ["group1", "group1", "group1"],
            name="Group",
            index=group_series.index,
        )
        filepath = str(clean_dir / "group_table.csv")
        group_series.to_csv(filepath, index=True)
        with pytest.raises(NotEnoughGroupsError):
            load_groups(filepath)


class TestGlyTraitInputData:
    @pytest.fixture
    def glycan_dict(self) -> dict[str, str]:
        return {"G1": "glycoct1", "G2": "glycoct2", "G3": "glycoct3"}

    def test_input_data(self, abundance_table, glycan_dict, group_series):
        data = GlyTraitInputData(
            abundance_table=abundance_table,
            glycans=glycan_dict,
            groups=group_series,
        )
        assert data.abundance_table.equals(abundance_table)
        assert data.glycans == {"G1": "glycoct1", "G2": "glycoct2", "G3": "glycoct3"}
        assert data.groups.equals(group_series)

    def test_no_series(self, abundance_table, glycan_dict):
        data = GlyTraitInputData(
            abundance_table=abundance_table,
            glycans=glycan_dict,
        )
        assert data.groups is None

    def test_glycans_not_same_as_in_abundance_table(self, abundance_table, glycan_dict):
        glycan_dict["G4"] = "glycoct4"
        with pytest.raises(FileFormatError):
            GlyTraitInputData(
                abundance_table=abundance_table,
                glycans=glycan_dict,
            )

    def test_samples_not_same_as_in_abundance_table(
        self, abundance_table, glycan_dict, group_series
    ):
        abundance_table.columns = ["sample1", "sample1", "sample3"]
        with pytest.raises(FileFormatError):
            GlyTraitInputData(
                abundance_table=abundance_table,
                glycans=glycan_dict,
                groups=group_series,
            )


@pytest.mark.usefixtures("patch_structure_from_string")
def test_load_input_data(
    abundance_table_file,
    group_series_file,
    structure_file,
    abundance_table,
    group_series,
):
    data = load_input_data(
        abundance_file=abundance_table_file,
        glycan_file=structure_file,
        mode="structure",
        group_file=group_series_file,
    )
    assert data.abundance_table.equals(abundance_table)
    assert data.groups.equals(group_series)
    assert data.glycans == {
        row["GlycanID"]: row["Structure"]
        for row in pd.read_csv(structure_file).to_dict(orient="records")
    }
