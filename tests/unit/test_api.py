import pytest
import pandas as pd

import glycoct as ct

from glytrait.api import Experiment, MissingDataError, InvalidOperationOrderError
from glytrait.exception import DataInputError


@pytest.fixture
def mock_auto_test(mocker):
    """Mock `glytrait.api.auto_test`."""
    mocker.patch("glytrait.api.auto_test", return_value=pd.DataFrame())


class TestExperiment:

    @pytest.fixture
    def abundance_file(self, clean_dir):
        """A temporary file of an abundance table."""
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
                "G1": [1.0, 1.0, 1.0],
                "G2": [1.0, 0.0, 0.0],
                "G3": [0.0, 0.0, 1.0],
            }
        )
        file = str(clean_dir / "abundance.csv")
        df.to_csv(file, index=False)
        return file

    @pytest.fixture
    def glycan_file(self, clean_dir):
        """A temporary file of a glycan table."""
        df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3"],
                "Structure": [ct.test_glycoct_1, ct.test_glycoct_2, ct.test_glycoct_3],
            }
        )
        file = str(clean_dir / "glycans.csv")
        df.to_csv(file, index=False)
        return file

    @pytest.fixture
    def mp_file(self, clean_dir):
        """A temporary file of a meta-property table."""
        df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3"],
                "nF": [1, 1, 0],
                "nS": [2, 0, 1],
            }
        )
        file = str(clean_dir / "meta_properties.csv")
        df.to_csv(file, index=False)
        return file

    @pytest.fixture
    def group_file(self, clean_dir):
        """A temporary file of a group table."""
        df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3"],
                "Group": ["A", "A", "B"],
            }
        )
        file = str(clean_dir / "groups.csv")
        df.to_csv(file, index=False)
        return file

    def test_init_with_abundance_and_glycan_files(self, abundance_file, glycan_file):
        """Test initiation with an abundance file and a glycan file."""
        exp = Experiment(abundance_file, glycan_file=glycan_file)

        assert exp.abundance_table is not None
        assert exp.meta_property_table is not None

    def test_init_with_abundance_and_mp_files(self, abundance_file, mp_file):
        """Test initiation with an abundance file and a meta-property file."""
        exp = Experiment(abundance_file, meta_property_file=mp_file)

        assert exp.abundance_table is not None
        assert exp.meta_property_table is not None

    def test_init_with_neither_glycan_or_mp_file(self, abundance_file):
        """Test initiation with neither a glycan file or a meta-property file."""
        with pytest.raises(DataInputError) as excinfo:
            Experiment(abundance_file)
        msg = (
            "At least one of `glycan_file` and `meta_property_file` should be provided."
        )
        assert str(excinfo.value) == msg

    def test_init_with_both_glycan_and_mp_file(
        self, abundance_file, glycan_file, mp_file
    ):
        """Test initiation with both a glycan file and a meta-property file."""
        with pytest.raises(DataInputError) as excinfo:
            Experiment(
                abundance_file, glycan_file=glycan_file, meta_property_file=mp_file
            )
        msg = "Only one of `glycan_file` and `meta_property_file` should be provided."
        assert str(excinfo.value) == msg

    def test_init_with_abundance_glycan_and_group_files(
        self, abundance_file, glycan_file, group_file
    ):
        """Test initiation with an abundance file, a glycan file, and a group file."""
        exp = Experiment(abundance_file, glycan_file=glycan_file, group_file=group_file)

        assert exp.abundance_table is not None
        assert exp.meta_property_table is not None
        assert exp.groups is not None

    def test_init_missing_glycans_in_glycan_file(self, abundance_file, clean_dir):
        """Test initiation when some glycans are missing in the glycan file."""
        glycan_df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G3"],  # G2 is missing
                "Structure": [ct.test_glycoct_1, ct.test_glycoct_3],
            }
        )
        glycan_file = str(clean_dir / "glycans.csv")
        glycan_df.to_csv(glycan_file, index=False)

        with pytest.raises(DataInputError) as excinfo:
            Experiment(abundance_file, glycan_file=glycan_file)
        msg = (
            "The following glycans in the abundance table do not have structures or "
            "compositions: G2."
        )
        assert str(excinfo.value) == msg

    def test_init_extra_glycans_in_glycan_file(self, abundance_file, clean_dir):
        """Test initiation when extra glycans exist in the glycan file."""
        glycan_df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3", "G4"],  # G4 is extra
                "Structure": [
                    ct.test_glycoct_1,
                    ct.test_glycoct_2,
                    ct.test_glycoct_3,
                    ct.test_glycoct_4,
                ],
            }
        )
        glycan_file = str(clean_dir / "glycans.csv")
        glycan_df.to_csv(glycan_file, index=False)

        # Extra glycans are ignored and should not raise an error
        exp = Experiment(abundance_file, glycan_file=glycan_file)

        assert exp.abundance_table is not None
        assert exp.meta_property_table is not None

    def test_init_with_missing_glycans_in_mp_file(self, abundance_file, clean_dir):
        """Test initiation when some glycans are missing in the meta-property file."""
        mp_df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G3"],  # missing G2
                "nF": [1, 0],
                "nS": [0, 1],
            }
        )
        mp_file = str(clean_dir / "meta_properties.csv")
        mp_df.to_csv(mp_file, index=False)

        with pytest.raises(DataInputError) as excinfo:
            Experiment(abundance_file, meta_property_file=mp_file)
        msg = "The following glycans in the abundance table do not have meta properties: G2."
        assert str(excinfo.value) == msg

    def test_init_with_extra_glycans_in_mp_file(self, abundance_file, clean_dir):
        """Test initiation when extra glycans exist in the meta-property file."""
        mp_df = pd.DataFrame(
            {
                "GlycanID": ["G1", "G2", "G3", "G4"],  # G4 is extra
                "nF": [1, 0, 1, 0],
                "nS": [0, 1, 0, 1],
            }
        )
        mp_file = str(clean_dir / "meta_properties.csv")
        mp_df.to_csv(mp_file, index=False)

        # Extra glycans are ignored and should not raise an error
        exp = Experiment(abundance_file, meta_property_file=mp_file)

        assert exp.abundance_table is not None
        assert exp.meta_property_table is not None

    def test_init_missing_samples_in_group_file(
        self, abundance_file, glycan_file, clean_dir
    ):
        """Test initiation when some samples are missing in the group file."""
        group_df = pd.DataFrame({"Sample": ["S1", "S3"], "Group": ["A", "B"]})
        group_file = str(clean_dir / "groups.csv")
        group_df.to_csv(group_file, index=False)

        with pytest.raises(DataInputError) as excinfo:
            Experiment(abundance_file, glycan_file=glycan_file, group_file=group_file)
        msg = "The following samples are in the abundance table but not in the groups: S2."
        assert str(excinfo.value) == msg

    def test_init_extra_samples_in_group_file(
        self, abundance_file, glycan_file, clean_dir
    ):
        """Test initiation when extra samples exist in the group file."""
        group_df = pd.DataFrame(
            {
                "Sample": ["S1", "S2", "S3", "S4"],
                "Group": ["A", "A", "B", "B"],
            }
        )
        group_file = str(clean_dir / "groups.csv")
        group_df.to_csv(group_file, index=False)

        with pytest.raises(DataInputError) as excinfo:
            Experiment(abundance_file, glycan_file=glycan_file, group_file=group_file)
        msg = "The following samples are in the groups but not in the abundance table: S4."
        assert str(excinfo.value) == msg

    @pytest.fixture
    def exp_no_groups(self, abundance_file, glycan_file):
        """An `Experiment` object without groups."""
        return Experiment(abundance_file, glycan_file=glycan_file)

    @pytest.fixture
    def exp_with_groups(self, abundance_file, glycan_file, group_file):
        """An `Experiment` object with groups."""
        return Experiment(
            abundance_file, glycan_file=glycan_file, group_file=group_file
        )

    def test_run_workflow_without_groups(self, exp_no_groups):
        """Test `run_workflow` when groups are not provided."""
        exp_no_groups.run_workflow()

        assert exp_no_groups.processed_abundance_table is not None
        assert exp_no_groups.derived_trait_table is not None
        assert exp_no_groups.filtered_derived_trait_table is not None
        with pytest.raises(MissingDataError):
            exp_no_groups.diff_results

    @pytest.mark.usefixtures("mock_auto_test")
    def test_run_workflow_with_groups(self, exp_with_groups):
        """Test `run_workflow` with groups provided."""
        exp_with_groups.run_workflow()

        assert exp_with_groups.processed_abundance_table is not None
        assert exp_with_groups.derived_trait_table is not None
        assert exp_with_groups.filtered_derived_trait_table is not None
        assert exp_with_groups.diff_results is not None

    def test_preprocess(self, exp_no_groups):
        """Test mannually calling `preprocess`."""
        exp_no_groups.preprocess()

        assert exp_no_groups.processed_abundance_table is not None
        with pytest.raises(MissingDataError):
            exp_no_groups.derived_trait_table

    def test_derive_traits(self, exp_no_groups):
        """Test mannually calling `derive_traits`."""
        exp_no_groups.preprocess()

        exp_no_groups.derive_traits()

        assert exp_no_groups.derived_trait_table is not None
        with pytest.raises(MissingDataError):
            exp_no_groups.filtered_derived_trait_table

    def test_post_filter(self, exp_no_groups):
        """Test mannually calling `post_filter`."""
        exp_no_groups.preprocess()
        exp_no_groups.derive_traits()

        exp_no_groups.post_filter()

        assert exp_no_groups.filtered_derived_trait_table is not None

    @pytest.mark.usefixtures("mock_auto_test")
    def test_diff_analysis(self, exp_with_groups):
        """Test mannually calling `diff_analysis`."""
        exp_with_groups.preprocess()
        exp_with_groups.derive_traits()
        exp_with_groups.post_filter()

        exp_with_groups.diff_analysis()

        assert exp_with_groups.diff_results is not None

    def test_derive_traits_before_preprocess(self, exp_no_groups):
        """Test calling `derived_traits` before `preprocess`."""
        with pytest.raises(InvalidOperationOrderError) as excinfo:
            exp_no_groups.derive_traits()
        assert str(excinfo.value) == "Please call the `preprocess` method first."

    def test_post_filter_before_derive_traits(self, exp_no_groups):
        """Test calling `post_filter` before `derive_traits`."""
        exp_no_groups.preprocess()

        with pytest.raises(InvalidOperationOrderError) as excinfo:
            exp_no_groups.post_filter()
        assert str(excinfo.value) == "Please call the `derive_traits` method first."

    @pytest.mark.usefixtures("mock_auto_test")
    def test_diff_analysis_before_post_filter(self, exp_with_groups):
        """Test calling `diff_analysis` before `post_filter`."""
        exp_with_groups.preprocess()
        exp_with_groups.derive_traits()

        with pytest.raises(InvalidOperationOrderError) as excinfo:
            exp_with_groups.diff_analysis()
        assert str(excinfo.value) == "Please call the `post_filter` method first."

    @pytest.mark.usefixtures("mock_auto_test")
    def test_diff_analysis_without_groups(self, exp_no_groups):
        """Test calling `diff_analysis` when groups are not provided."""
        exp_no_groups.preprocess()
        exp_no_groups.derive_traits()
        exp_no_groups.post_filter()

        with pytest.raises(MissingDataError) as excinfo:
            exp_no_groups.diff_analysis()
        assert (
            str(excinfo.value)
            == "Group information is required for differential analysis."
        )

    def test_get_diff_results_without_groups(self, exp_no_groups):
        """Test getting `diff_results` when groups are not provided."""
        exp_no_groups.run_workflow()

        with pytest.raises(MissingDataError) as excinfo:
            exp_no_groups.diff_results
        msg = (
            "Grouping information is missing. Please create a new `Experiment` "
            "object with `group_file` first."
        )
        assert str(excinfo.value) == msg

    def test_modify_abundance(self, exp_no_groups):
        """Test if modifying the abundance table doesn't work."""
        exp_no_groups.abundance_table.drop(columns=["G1"], inplace=True)

        assert exp_no_groups.abundance_table.shape == (3, 3)

    def test_modify_groups(self, exp_with_groups):
        """Test if modifying the group table doesn't work."""
        exp_with_groups.groups.drop(labels=["S1"], inplace=True)

        assert exp_with_groups.groups.shape == (3,)

    def test_modify_meta_properties(self, exp_no_groups):
        """Test if modifying the meta property table doesn't work."""
        exp_no_groups.meta_property_table.drop(index=["G1"], inplace=True)

        assert exp_no_groups.meta_property_table.shape[0] == 3

    def test_modify_processed_abundance(self, exp_no_groups):
        """Test if modifying the processed abundance table doesn't work."""
        exp_no_groups.preprocess()
        exp_no_groups.processed_abundance_table.drop(columns=["G1"], inplace=True)

        assert exp_no_groups.processed_abundance_table.shape == (3, 3)

    def test_modify_derived_traits(self, exp_no_groups):
        """Test if modifying the derived trait table doesn't work."""
        exp_no_groups.preprocess()
        exp_no_groups.derive_traits()
        exp_no_groups.derived_trait_table.drop(index=["S1"], inplace=True)

        assert exp_no_groups.derived_trait_table.shape[0] == 3

    def test_modify_filtered_derived_traits(self, exp_no_groups):
        """Test if modifying the filtered derived trait table doesn't work."""
        exp_no_groups.preprocess()
        exp_no_groups.derive_traits()
        exp_no_groups.post_filter()
        exp_no_groups.filtered_derived_trait_table.drop(index=["S1"], inplace=True)

        assert exp_no_groups.filtered_derived_trait_table.shape[0] == 3

    def test_modify_diff_results(self, exp_with_groups):
        """Test if modifying the differential analysis results doesn't work."""
        exp_with_groups.run_workflow()
        exp_with_groups.diff_results["new"] = pd.DataFrame()

        assert "new_col" not in exp_with_groups.diff_results

    @pytest.mark.usefixtures("mock_auto_test")
    @pytest.mark.parametrize(
        "method",
        [
            "preprocess",
            "derive_traits",
            "post_filter",
            "diff_analysis",
        ],
    )
    def test_rerunning_step_methods_clears_afterward_results(
        self, exp_with_groups, method
    ):
        """Test if re-running a step method clears the afterward results."""
        method_names = ["preprocess", "derive_traits", "post_filter", "diff_analysis"]
        methods = [
            exp_with_groups.preprocess,
            exp_with_groups.derive_traits,
            exp_with_groups.post_filter,
            exp_with_groups.diff_analysis,
        ]
        attributes = {
            "preprocess": "processed_abundance_table",
            "derive_traits": "derived_trait_table",
            "post_filter": "filtered_derived_trait_table",
            "diff_analysis": "diff_results",
        }
        method_index = method_names.index(method)
        exp_with_groups.run_workflow()

        methods[method_index]()

        for afterward_method_index in range(method_index + 1, 4):
            method_name = method_names[afterward_method_index]
            with pytest.raises(MissingDataError):
                getattr(exp_with_groups, attributes[method_name])

    def test_try_multiple_formulas(self, exp_no_groups):
        """Test `try_formula` with multiple formula expressions."""
        exp_no_groups.preprocess()
        formulas = [
            "TC = [type == 'complex'] / [1]",
            "TS = [nS > 0] / [1]",
        ]

        result = exp_no_groups.try_formulas(formulas)

        assert result.shape == (3, 2)

    def test_try_one_formula_no_squeeze(self, exp_no_groups):
        """Test `try_formula` with one formula expression with `squeeze` setting to False."""
        exp_no_groups.preprocess()

        result = exp_no_groups.try_formulas("TS = [nS > 0] / [1]", squeeze=False)

        assert isinstance(result, pd.DataFrame)

    def test_try_one_formula_squeeze(self, exp_no_groups):
        """Test `try_formula` with one formula expression with `squeeze` setting to True."""
        exp_no_groups.preprocess()

        resource = exp_no_groups.try_formulas("TS = [nS > 0] / [1]", squeeze=True)

        assert isinstance(resource, pd.Series)

    def test_try_one_formula_default_squeeze(self, exp_no_groups):
        """Test `try_formula` with one formula expression with `squeeze` setting to True."""
        exp_no_groups.preprocess()

        resource = exp_no_groups.try_formulas("TS = [nS > 0] / [1]")

        assert isinstance(resource, pd.Series)
