import pandas as pd
import pytest
from numpy import dtype

import glytrait
from glytrait.exception import (
    DataInputError,
)
from glytrait.glycan import StructureDict
from glytrait.load_data import (
    load_input_data,
    AbundanceCSVLoader,
    GroupsCSVLoader,
    GlycanCSVLoader,
    DFValidator,
)


class TestReadCSV:
    """Test `_read_csv`."""

    def test_basic(self, clean_dir):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        filepath = clean_dir / "test.csv"
        df.to_csv(filepath, index=False)
        result = glytrait.load_data._read_csv(filepath)
        pd.testing.assert_frame_equal(result, df)

    def test_file_not_found(self, clean_dir):
        filepath = clean_dir / "abundance_table.csv"
        loader = AbundanceCSVLoader(filepath=filepath)
        with pytest.raises(FileNotFoundError):
            loader.load()

    def test_wrong_format(self, clean_dir):
        content = "G1,G2,G3\n1,2,3\n4,5,6\n7,8,9,10"
        filepath = clean_dir / "abundance_table.csv"
        filepath.write_text(content)
        loader = AbundanceCSVLoader(filepath=filepath)
        with pytest.raises(DataInputError) as excinfo:
            loader.load()
        assert "This CSV file could not be parsed." in str(excinfo.value)

    def test_empty_file(self, clean_dir):
        filepath = clean_dir / "abundance_table.csv"
        filepath.touch()
        loader = AbundanceCSVLoader(filepath=filepath)
        with pytest.raises(DataInputError) as excinfo:
            loader.load()
        assert f"Empty CSV file." in str(excinfo.value)


class TestGroupsCSVLoader:
    def test_basic(self, mocker):
        df = pd.DataFrame(
            {
                "Sample": ["sample1", "sample2", "sample3"],
                "Group": ["group1", "group2", "group3"],
            },
        )
        mocker.patch("glytrait.load_data._read_csv", return_value=df)
        loader = GroupsCSVLoader(filepath="group_table.csv", validator=lambda df: None)
        result = loader.load()
        expected = pd.Series(
            ["group1", "group2", "group3"],
            name="Group",
            index=pd.Index(["sample1", "sample2", "sample3"], name="Sample"),
        )
        pd.testing.assert_series_equal(result, expected)
        glytrait.load_data._read_csv.assert_called_once_with("group_table.csv")

    def test_default_validator(self):
        loader = GroupsCSVLoader(filepath="group_table.csv")
        assert loader.validator == DFValidator(
            must_have=["Group", "Sample"], unique=["Sample"], types={"Sample": "object"}
        )


class TestAbundanceCSVLoader:
    def test_basic(self, mocker):
        df = pd.DataFrame(
            {
                "Sample": ["sample1", "sample2", "sample3"],
                "G1": [1, 2, 3],
                "G2": [4, 5, 6],
            }
        )
        mocker.patch("glytrait.load_data._read_csv", return_value=df)
        loader = AbundanceCSVLoader(
            filepath="abundance_table.csv", validator=lambda df: None
        )
        result = loader.load()
        expected = pd.DataFrame(
            {
                "G1": [1, 2, 3],
                "G2": [4, 5, 6],
            },
            index=pd.Index(["sample1", "sample2", "sample3"], name="Sample"),
        )
        pd.testing.assert_frame_equal(result, expected)
        glytrait.load_data._read_csv.assert_called_once_with("abundance_table.csv")

    def test_default_validator(self):
        loader = AbundanceCSVLoader(filepath="abundance_table.csv")
        assert loader.validator == DFValidator(
            must_have=["Sample"], unique=["Sample"], types={"Sample": "object"}
        )


class TestGlycanCSVLoader:
    @staticmethod
    def fake_parser(it):
        return StructureDict({name: string for name, string in it})

    @pytest.mark.parametrize("mode", ["structure", "composition"])
    def test_basic(self, mocker, mode):
        glycan_col = "Structure" if mode == "structure" else "Composition"
        df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3"],
                glycan_col: ["glycan1", "glycan2", "glycan3"],
            },
        )
        mocker.patch("glytrait.load_data._read_csv", return_value=df)
        loader = GlycanCSVLoader(
            "structure.csv",
            mode=mode,
            parser=self.fake_parser,
            validator=lambda df: None,
        )
        result = loader.load()
        expected = {"G1": "glycan1", "G2": "glycan2", "G3": "glycan3"}
        assert result == expected

    @pytest.mark.parametrize(
        "mode, parser",
        [
            ("structure", glytrait.glycan.parse_structures),
            ("composition", glytrait.glycan.parse_compositions),
        ],
    )
    def test_default_parser(self, mode, parser):
        loader = GlycanCSVLoader("structure.csv", mode=mode)
        assert loader.parser == parser

    @pytest.mark.parametrize(
        "mode, glycan_col",
        [
            ("structure", "Structure"),
            ("composition", "Composition"),
        ],
    )
    def test_default_validator(self, mode, glycan_col):
        loader = GlycanCSVLoader("structure.csv", mode=mode)
        assert loader.validator == DFValidator(
            must_have=["GlycanID", glycan_col],
            unique=["GlycanID", glycan_col],
            types={"GlycanID": "object", glycan_col: "object"},
        )


class TestDFValidator:
    def test_no_validation(self):
        validator = DFValidator()
        assert validator(pd.DataFrame({"A": [1, 2], "B": [3, 4]})) is None

    def test_has_required_columns(self):
        validator = DFValidator(must_have=["A", "B"])
        assert validator(pd.DataFrame({"A": [1, 2], "B": [3, 4]})) is None

    def test_missing_required_columns(self):
        validator = DFValidator(must_have=["A", "B"])
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [1, 2], "C": [3, 4]}))
        assert "The following columns are missing: B" in str(excinfo.value)

    def test_unique_columns(self):
        validator = DFValidator(unique=["A"])
        assert validator(pd.DataFrame({"A": [1, 2], "B": [3, 3]})) is None

    def test_specified_columns_not_unique(self):
        validator = DFValidator(unique=["A"])
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [1, 1], "B": [3, 3]}))
        assert "The following columns are not unique: A" in str(excinfo.value)

    def test_unique_columns_not_exist(self):
        validator = DFValidator(unique=["A"])
        assert validator(pd.DataFrame({"B": [3, 3]})) is None

    def test_type_check_pass(self):
        validator = DFValidator(types={"A": dtype("int64"), "B": "object"})
        assert validator(pd.DataFrame({"A": [1, 2], "B": ["3", "4"]})) is None

    def test_type_check_fail(self):
        validator = DFValidator(types={"A": dtype("int64"), "B": "object"})
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [1, 2], "B": [3, 4]}))
        msg = (
            "The following columns have incorrect types: B. "
            "Expected types: {'B': 'object'}, got: {'B': dtype('int64')}."
        )
        assert msg in str(excinfo.value)

    def test_type_check_missing_column(self):
        validator = DFValidator(types={"A": dtype("int64"), "B": "object"})
        assert validator(pd.DataFrame({"A": [1, 2]})) is None

    def test_default_type_check_pass(self):
        validator = DFValidator(default_type="object")
        assert validator(pd.DataFrame({"A": ["3", "4"], "B": ["3", "4"]})) is None

    def test_default_type_check_fail(self):
        validator = DFValidator(default_type=dtype("float64"))
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [3.1, 3.2], "B": [3, 4]}))
        msg = (
            "The following columns have incorrect types: B. "
            "Expected types: float64, got: {'B': dtype('int64')}."
        )
        assert msg in str(excinfo.value)
