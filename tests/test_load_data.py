import pandas as pd
import pytest

from glytrait.load_data import load_abundance_table, load_structures
from glytrait.exception import FileTypeError, FileFormatError, StructureParseError


@pytest.fixture
def abundance_table() -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "sample1": [1, 2, 3],
            "sample2": [4, 5, 6],
            "sample3": [7, 8, 9],
        },
        index=pd.Index(["G1", "G2", "G3"], name="GlycanID"),
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


@pytest.mark.usefixtures("patch_structure_from_string")
class TestLoadStructures:
    def test_basic(self, structure_file):
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
