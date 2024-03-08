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
    CSVLoader,
    AbundanceCSVLoader,
    GroupsCSVLoader,
    GlycanCSVLoader,
    DFValidator,
    GlyTraitInputData,
    InputDataValidator,
    load_input_data,
    load_input_data_from_csv,
)


class TestCSVLoader:

    def test_read_csv(self, clean_dir):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        filepath = clean_dir / "test.csv"
        df.to_csv(filepath, index=False)
        result = CSVLoader(filepath).read_csv()
        pd.testing.assert_frame_equal(result, df)

    def test_read_csv_file_not_found(self, clean_dir):
        filepath = clean_dir / "abundance_table.csv"
        with pytest.raises(FileNotFoundError):
            result = CSVLoader(filepath).read_csv()

    def test_read_csv_wrong_format(self, clean_dir):
        content = "G1,G2,G3\n1,2,3\n4,5,6\n7,8,9,10"
        filepath = clean_dir / "abundance_table.csv"
        filepath.write_text(content)
        with pytest.raises(DataInputError) as excinfo:
            result = CSVLoader(filepath).read_csv()
        assert "This CSV file could not be parsed." in str(excinfo.value)

    def test_read_csv_empty_file(self, clean_dir):
        filepath = clean_dir / "abundance_table.csv"
        filepath.touch()
        with pytest.raises(DataInputError) as excinfo:
            result = CSVLoader(filepath).read_csv()
        assert f"Empty CSV file." in str(excinfo.value)

    def test_load_df(self, mocker):
        validator = mocker.Mock()
        returned_df = mocker.Mock()
        loader = CSVLoader(
            filepath="test.csv", reader=lambda x: returned_df, validator=validator
        )
        result = loader.load_df()
        assert result == returned_df


class TestGroupsCSVLoader:
    def test_basic(self):
        df = pd.DataFrame(
            {
                "Sample": ["sample1", "sample2", "sample3"],
                "Group": ["group1", "group2", "group3"],
            },
        )
        loader = GroupsCSVLoader(
            filepath="group_table.csv", reader=lambda x: df, validator=lambda df: None
        )
        result = loader.load()
        expected = pd.Series(
            ["group1", "group2", "group3"],
            name="Group",
            index=pd.Index(["sample1", "sample2", "sample3"], name="Sample"),
        )
        pd.testing.assert_series_equal(result, expected)

    def test_default_validator(self):
        loader = GroupsCSVLoader(filepath="group_table.csv")
        assert loader.validator == DFValidator(
            must_have=["Group", "Sample"], unique=["Sample"], types={"Sample": "object"}
        )


class TestAbundanceCSVLoader:
    def test_basic(self):
        df = pd.DataFrame(
            {
                "Sample": ["sample1", "sample2", "sample3"],
                "G1": [1, 2, 3],
                "G2": [4, 5, 6],
            }
        )
        loader = AbundanceCSVLoader(
            filepath="abundance_table.csv",
            reader=lambda x: df,
            validator=lambda df: None,
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
    def test_basic(self, mode):
        glycan_col = "Structure" if mode == "structure" else "Composition"
        df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3"],
                glycan_col: ["glycan1", "glycan2", "glycan3"],
            },
        )
        loader = GlycanCSVLoader(
            "structure.csv",
            mode=mode,
            parser=self.fake_parser,
            reader=lambda x: df,
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


class TestInputDataValidator:

    def test_registered_validators(self):
        validators = set(InputDataValidator.validators)
        expected = {
            glytrait.load_data.same_samples_in_abundance_and_groups,
            glytrait.load_data.all_glycans_have_structures_or_compositions,
        }
        assert validators == expected


class TestSameSamplesInAbundanceAndGroups:

    @staticmethod
    def make_abund_df(samples: list[str]):
        """Helper function to create a simple abundance table with given sample names."""
        return pd.DataFrame(
            {"G1": [1] * len(samples), "G2": [1] * len(samples)}, index=samples
        )

    @staticmethod
    def make_group_s(samples: list[str]):
        """Helper function to create a simple group series with given sample names."""
        return pd.Series(
            ["group1"] * len(samples), index=pd.Index(samples, name="Sample")
        )

    def test_same(self):
        abundance = self.make_abund_df(["sample1", "sample2"])
        groups = self.make_group_s(["sample1", "sample2"])
        data = GlyTraitInputData(abundance_table=abundance, groups=groups, glycans=None)
        glytrait.load_data.same_samples_in_abundance_and_groups(data)

    def test_sample_missing_in_groups(self):
        abundance = self.make_abund_df(["sample1", "sample2"])
        groups = self.make_group_s(["sample1"])
        data = GlyTraitInputData(abundance_table=abundance, groups=groups, glycans=None)
        with pytest.raises(DataInputError) as excinfo:
            glytrait.load_data.same_samples_in_abundance_and_groups(data)
        msg = "The following samples are in the abundance table but not in the groups: sample2."
        assert msg in str(excinfo.value)

    def test_sample_missing_in_abundance(self):
        abundance = self.make_abund_df(["sample1"])
        groups = self.make_group_s(["sample1", "sample2"])
        data = GlyTraitInputData(abundance_table=abundance, groups=groups, glycans=None)
        with pytest.raises(DataInputError) as excinfo:
            glytrait.load_data.same_samples_in_abundance_and_groups(data)
        msg = "The following samples are in the groups but not in the abundance table: sample2."
        assert msg in str(excinfo.value)

    def test_sample_missing_in_both(self):
        abundance = self.make_abund_df(["sample1"])
        groups = self.make_group_s(["sample2"])
        data = GlyTraitInputData(abundance_table=abundance, groups=groups, glycans=None)
        with pytest.raises(DataInputError) as excinfo:
            glytrait.load_data.same_samples_in_abundance_and_groups(data)
        msg1 = "The following samples are in the abundance table but not in the groups: sample1."
        msg2 = "The following samples are in the groups but not in the abundance table: sample2."
        assert msg1 in str(excinfo.value)
        assert msg2 in str(excinfo.value)


class TestAllGlycansHaveStructuresOrComposition:

    @staticmethod
    def make_abund_df(glycans: list[str]) -> pd.DataFrame:
        """Helper function to create a simple abundance table with given glycan names."""
        return pd.DataFrame({glycan: [1] for glycan in glycans}, index=["sample1"])

    def test_all_have_structures(self):
        glycans = {"G1": "glycan1", "G2": "glycan2"}
        abund_df = self.make_abund_df(glycans.keys())
        data = GlyTraitInputData(abundance_table=abund_df, groups=None, glycans=glycans)
        glytrait.load_data.all_glycans_have_structures_or_compositions(data)

    def test_missing_structures(self):
        glycans = {"G1": "glycan1"}
        abund_df = self.make_abund_df(["G1", "G2"])
        data = GlyTraitInputData(abundance_table=abund_df, groups=None, glycans=glycans)
        with pytest.raises(DataInputError) as excinfo:
            glytrait.load_data.all_glycans_have_structures_or_compositions(data)
        msg = (
            f"The following glycans in the abundance table do not have structures or "
            f"compositions: G2."
        )
        assert msg in str(excinfo.value)

    def test_glycans_in_structures_not_in_abundance(self):
        glycans = {"G1": "glycan1", "G2": "glycan2"}
        abund_df = self.make_abund_df(["G1"])
        data = GlyTraitInputData(abundance_table=abund_df, groups=None, glycans=glycans)
        assert (
            glytrait.load_data.all_glycans_have_structures_or_compositions(data) is None
        )


def test_load_input_data(mocker):
    abundance_table = mocker.Mock()
    glycan_dict = mocker.Mock()
    group_seris = mocker.Mock()

    abundance_loader = mocker.Mock()
    abundance_loader.load.return_value = abundance_table
    glycan_loader = mocker.Mock()
    glycan_loader.load.return_value = glycan_dict
    group_loader = mocker.Mock()
    group_loader.load.return_value = group_seris

    validator = mocker.Mock()

    result = load_input_data(
        abundance_loader=abundance_loader,
        glycan_loader=glycan_loader,
        group_loader=group_loader,
        validator=validator,
    )

    assert result.abundance_table == abundance_table
    assert result.glycans == glycan_dict
    assert result.groups == group_seris

    abundance_loader.load.assert_called_once()
    glycan_loader.load.assert_called_once()
    group_loader.load.assert_called_once()
    validator.assert_called_once_with(result)


@pytest.mark.parametrize("mode", ["structure", "composition"])
@pytest.mark.parametrize("group_file", [None, "group.csv"])
def test_load_input_data_from_csv(mocker, mode, group_file):
    input_data = mocker.Mock()
    mocker.patch(
        "glytrait.load_data.load_input_data", autospec=True, return_value=input_data
    )
    abund_loader_mock = mocker.patch(
        "glytrait.load_data.AbundanceCSVLoader", autospec=True
    )
    glycan_loader_mock = mocker.patch(
        "glytrait.load_data.GlycanCSVLoader", autospec=True
    )
    group_loader_mock = mocker.patch(
        "glytrait.load_data.GroupsCSVLoader", autospec=True
    )

    result = load_input_data_from_csv(
        abundance_file="abundance.csv",
        glycan_file="glycan.csv",
        group_file=group_file,
        mode=mode,
    )

    assert result == input_data
    abund_loader_mock.assert_called_once_with(filepath="abundance.csv")
    glycan_loader_mock.assert_called_once_with(filepath="glycan.csv", mode=mode)
    if group_file:
        group_loader_mock.assert_called_once_with(filepath="group.csv")
    else:
        group_loader_mock.assert_not_called()

    params = {
        "abundance_loader": abund_loader_mock.return_value,
        "glycan_loader": glycan_loader_mock.return_value,
        "group_loader": group_loader_mock.return_value,
    }
    if not group_file:
        params["group_loader"] = None
    glytrait.load_data.load_input_data.assert_called_once_with(**params)
