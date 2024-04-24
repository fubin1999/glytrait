import pytest
import pandas as pd

from glytrait.data_input import load_data
from glytrait.api import Experiment


class TestExperiment:

    @pytest.fixture
    def input_data(self):
        abund_df = pd.read_csv("tests/integration/data/abundance.csv")
        glycan_df = pd.read_csv("tests/integration/data/structures.csv")
        group_df = pd.read_csv("tests/integration/data/groups.csv")
        return load_data(abund_df, glycan_df, group_df, mode="structure")

    @pytest.fixture
    def experiment(self, input_data):
        return Experiment(input_data)

    def test_whole_workflow(self, experiment):
        experiment.preprocess()
        experiment.extract_meta_properties()
        experiment.derive_traits()
        experiment.post_filter()
        experiment.diff_analysis()

        assert isinstance(experiment.processed_abundance_table, pd.DataFrame)
        assert isinstance(experiment.meta_property_table, pd.DataFrame)
        assert isinstance(experiment.derived_trait_table, pd.DataFrame)
        assert isinstance(experiment.filtered_derived_trait_table, pd.DataFrame)
        assert isinstance(experiment.diff_results, dict)
