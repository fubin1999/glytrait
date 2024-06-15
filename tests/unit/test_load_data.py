import pandas as pd
import pytest
from numpy import dtype

from glytrait.exception import DataInputError
import glytrait.data_input as di


class TestDFValidator:
    def test_no_validation(self):
        validator = di.DFValidator()
        assert validator(pd.DataFrame({"A": [1, 2], "B": [3, 4]})) is None

    def test_has_required_columns(self):
        validator = di.DFValidator(must_have=["A", "B"])
        assert validator(pd.DataFrame({"A": [1, 2], "B": [3, 4]})) is None

    def test_missing_required_columns(self):
        validator = di.DFValidator(must_have=["A", "B"])
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [1, 2], "C": [3, 4]}))
        assert "The following columns are missing: B" in str(excinfo.value)

    def test_unique_columns(self):
        validator = di.DFValidator(unique=["A"])
        assert validator(pd.DataFrame({"A": [1, 2], "B": [3, 3]})) is None

    def test_specified_columns_not_unique(self):
        validator = di.DFValidator(unique=["A"])
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [1, 1], "B": [3, 3]}))
        assert "The following columns are not unique: A" in str(excinfo.value)

    def test_unique_columns_not_exist(self):
        validator = di.DFValidator(unique=["A"])
        assert validator(pd.DataFrame({"B": [3, 3]})) is None

    def test_type_check_pass(self):
        validator = di.DFValidator(types={"A": dtype("int64"), "B": "object"})
        assert validator(pd.DataFrame({"A": [1, 2], "B": ["3", "4"]})) is None

    def test_type_check_fail(self):
        validator = di.DFValidator(types={"A": dtype("int64"), "B": "object"})
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [1, 2], "B": [3, 4]}))
        msg = (
            "The following columns have incorrect types: B. "
            "Expected types: {'B': 'object'}, got: {'B': dtype('int64')}."
        )
        assert msg in str(excinfo.value)

    def test_type_check_missing_column(self):
        validator = di.DFValidator(types={"A": dtype("int64"), "B": "object"})
        assert validator(pd.DataFrame({"A": [1, 2]})) is None

    def test_default_type_check_pass(self):
        validator = di.DFValidator(default_type="object")
        assert validator(pd.DataFrame({"A": ["3", "4"], "B": ["3", "4"]})) is None

    def test_default_type_check_fail(self):
        validator = di.DFValidator(default_type=dtype("float64"))
        with pytest.raises(DataInputError) as excinfo:
            validator(pd.DataFrame({"A": [3.1, 3.2], "B": [3, 4]}))
        msg = (
            "The following columns have incorrect types: B. "
            "Expected types: float64, got: {'B': dtype('int64')}."
        )
        assert msg in str(excinfo.value)


class TestLoadAbundance:

    def test_basic(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
                "G1": [1.0, 2.0, 3.0],
                "G2": [4.0, 5.0, 6.0],
                "G3": [7.0, 8.0, 9.0],
            }
        )
        result = di.load_abundance(df)
        expected = pd.DataFrame(
            {
                "G1": [1.0, 2.0, 3.0],
                "G2": [4.0, 5.0, 6.0],
                "G3": [7.0, 8.0, 9.0],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_no_sample_col(self):
        df = pd.DataFrame(
            {
                "G1": [1.0, 2.0, 3.0],
                "G2": [4.0, 5.0, 6.0],
                "G3": [7.0, 8.0, 9.0],
            }
        )
        with pytest.raises(DataInputError):
            di.load_abundance(df)

    def test_duplicate_samples(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S1", "S3"],
                "G1": [1.0, 2.0, 3.0],
                "G2": [4.0, 5.0, 6.0],
                "G3": [7.0, 8.0, 9.0],
            }
        )
        with pytest.raises(DataInputError):
            di.load_abundance(df)

    def test_wrong_type(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
                "G1": ["A", "B", "C"],
                "G2": [4.0, 5.0, 6.0],
                "G3": [7.0, 8.0, 9.0],
            }
        )
        with pytest.raises(DataInputError):
            di.load_abundance(df)

    @pytest.mark.skip("Not implemented")
    def test_only_sample_col(self):
        df = pd.DataFrame({"Sample": ["S1", "S2", "S3"]})
        with pytest.raises(DataInputError):
            di.load_abundance(df)


class TestLoadGroups:

    def test_basic(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
                "Group": ["G1", "G2", "G3"],
            }
        )
        result = di.load_groups(df)
        expected = pd.Series(
            ["G1", "G2", "G3"],
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
            name="Group",
        )
        pd.testing.assert_series_equal(result, expected)

    @pytest.mark.skip("Not implemented")
    def test_wrong_cols(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
                "Group": ["G1", "G2", "G3"],
                "Extra": [1, 2, 3],
            }
        )
        with pytest.raises(DataInputError):
            di.load_groups(df)

    def test_missing_sample_col(self):
        df = pd.DataFrame(
            {
                "Group": ["G1", "G2", "G3"],
            }
        )
        with pytest.raises(DataInputError):
            di.load_groups(df)

    def test_missing_group_col(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
            }
        )
        with pytest.raises(DataInputError):
            di.load_groups(df)

    def test_duplicate_samples(self):
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S1", "S3"],
                "Group": ["G1", "G2", "G3"],
            }
        )
        with pytest.raises(DataInputError):
            di.load_groups(df)


class TestLoadGlycans:

    def test_basic(self):
        df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3"],
                "Structure": ["A", "B", "C"],
            }
        )
        result = di.load_glycans(
            df, mode="structure", parser=lambda x: {k: v for k, v in x}
        )
        expected = {"G1": "A", "G2": "B", "G3": "C"}
        assert result == expected


class TestLoadMetaProperties:

    def test_basic(self):
        df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3"],
                "MP1": [1, 2, 3],
                "MP2": ["A", "B", "A"],
            }
        )
        result = di.load_meta_property(df)
        expected = pd.DataFrame(
            {
                "MP1": [1, 2, 3],
                "MP2": ["A", "B", "A"],
            },
            index=pd.Index(["G1", "G2", "G3"], name="GlycanID")
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_missing_glycan_id_col(self):
        df = pd.DataFrame(
            {
                "MP1": [1, 2, 3],
                "MP2": ["A", "B", "A"],
            }
        )
        with pytest.raises(DataInputError):
            di.load_meta_property(df)


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
        di.check_same_samples_in_abund_and_groups(abundance, groups)

    def test_sample_missing_in_groups(self):
        abundance = self.make_abund_df(["sample1", "sample2"])
        groups = self.make_group_s(["sample1"])
        with pytest.raises(DataInputError) as excinfo:
            di.check_same_samples_in_abund_and_groups(abundance, groups)
        msg = "The following samples are in the abundance table but not in the groups: sample2."
        assert msg in str(excinfo.value)

    def test_sample_missing_in_abundance(self):
        abundance = self.make_abund_df(["sample1"])
        groups = self.make_group_s(["sample1", "sample2"])
        with pytest.raises(DataInputError) as excinfo:
            di.check_same_samples_in_abund_and_groups(abundance, groups)
        msg = "The following samples are in the groups but not in the abundance table: sample2."
        assert msg in str(excinfo.value)

    def test_sample_missing_in_both(self):
        abundance = self.make_abund_df(["sample1"])
        groups = self.make_group_s(["sample2"])
        with pytest.raises(DataInputError) as excinfo:
            di.check_same_samples_in_abund_and_groups(abundance, groups)
        msg1 = "The following samples are in the abundance table but not in the groups: sample1."
        msg2 = "The following samples are in the groups but not in the abundance table: sample2."
        assert msg1 in str(excinfo.value)
        assert msg2 in str(excinfo.value)


def make_abund_df(glycans: list[str]) -> pd.DataFrame:
    """Helper function to create a simple abundance table with given glycan names."""
    return pd.DataFrame({glycan: [1] for glycan in glycans}, index=["sample1"])


class TestAllGlycansHaveStructuresOrComposition:

    def test_all_have_structures(self):
        glycans = {"G1": "glycan1", "G2": "glycan2"}
        abund_df = make_abund_df(glycans.keys())
        di.check_all_glycans_have_struct_or_comp(abund_df, glycans)

    def test_missing_structures(self):
        glycans = {"G1": "glycan1"}
        abund_df = make_abund_df(["G1", "G2"])
        with pytest.raises(DataInputError) as excinfo:
            di.check_all_glycans_have_struct_or_comp(abund_df, glycans)
        msg = (
            "The following glycans in the abundance table do not have structures or "
            "compositions: G2."
        )
        assert msg in str(excinfo.value)

    def test_glycans_in_structures_not_in_abundance(self):
        glycans = {"G1": "glycan1", "G2": "glycan2"}
        abund_df = make_abund_df(["G1"])
        di.check_all_glycans_have_struct_or_comp(abund_df, glycans)


def make_mp_df(glycans: list[str]) -> pd.DataFrame:
    """Helper function to create a simple mp table with the given glycans names."""
    return pd.DataFrame(
        {
            "MP1": [1] * len(glycans),
            "MP2": ["a"] * len(glycans)
        },
        index=pd.Index(glycans, name="GlycanID")
    )


class TestAllGlycansHaveMP:

    def test_all_have_structures(self):
        abund_df = make_abund_df(["G1", "G2"])
        mp_df = make_mp_df(["G1", "G2"])
        di.check_all_glycans_have_mp(abund_df, mp_df)

    def test_missing_structures(self):
        abund_df = make_abund_df(["G1", "G2"])
        mp_df = make_mp_df(["G1"])
        with pytest.raises(DataInputError) as excinfo:
            di.check_all_glycans_have_mp(abund_df, mp_df)
        msg = "The following glycans in the abundance table do not have meta properties: G2."
        assert msg in str(excinfo.value)

    def test_glycans_in_structures_not_in_abundance(self):
        abund_df = make_abund_df(["G1"])
        mp_df = make_mp_df(["G1", "G2"])
        di.check_all_glycans_have_mp(abund_df, mp_df)


class TestGlyTraitInputData:

    @pytest.fixture(autouse=True)
    def patch_checker(self, mocker):
        mocker.patch("glytrait.data_input.check_same_samples_in_abund_and_groups")
        mocker.patch("glytrait.data_input.check_all_glycans_have_struct_or_comp")

    @pytest.fixture
    def abundance_table(self):
        return pd.DataFrame(
            {
                "G1": [1.0, 2.0, 3.0],
                "G2": [4.0, 5.0, 6.0],
            },
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
        )

    @pytest.fixture
    def glycans(self):
        return {
            "G1": "glycan1",
            "G2": "glycan2",
        }

    @pytest.fixture
    def groups(self):
        return pd.Series(
            ["A", "B", "A"],
            index=pd.Index(["S1", "S2", "S3"], name="Sample"),
            name="Group",
        )

    @pytest.fixture
    def input_data(self, abundance_table, glycans, groups):
        return di.GlyTraitInputData(
            abundance_table=abundance_table,
            glycans=glycans,
            groups=groups,
        )

    def test_abundance_table_getter(self, input_data, abundance_table):
        # Test if the getter returns the correct value
        assert input_data.abundance_table.equals(abundance_table)

        # Test if a copy is returned
        input_data.abundance_table["New"] = 1  # Try to modify the returned value
        assert input_data.abundance_table.equals(
            abundance_table
        )  # Check if the original value is not modified

    def test_abundance_table_setter(self, input_data):
        input_data.abundance_table = pd.DataFrame()
        assert input_data.abundance_table.empty
        assert di.check_same_samples_in_abund_and_groups.called_once()
        assert di.check_all_glycans_have_struct_or_comp.called_once()

    def test_glycans_getter(self, input_data, glycans):
        # Test if the getter returns the correct value
        assert input_data.glycans == glycans

        # Test if a copy is returned
        input_data.glycans["G3"] = "glycan3"
        assert input_data.glycans == glycans

    def test_glycans_setter(self, input_data):
        input_data.glycans = {}
        assert input_data.glycans == {}
        assert di.check_all_glycans_have_struct_or_comp.called_once()

    def test_groups_getter(self, input_data, groups):
        # Test if the getter returns the correct value
        assert input_data.groups.equals(groups)

        # Test if a copy is returned
        input_data.groups["S4"] = "C"
        assert input_data.groups.equals(groups)

    def test_groups_setter(self, input_data):
        input_data.groups = pd.Series()
        assert input_data.groups.empty
        assert di.check_same_samples_in_abund_and_groups.called_once()
