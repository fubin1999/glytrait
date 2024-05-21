import pytest
import pandas as pd

from glytrait.api import Experiment


class TestExperiment:

    @pytest.fixture
    def experiment(self):
        return Experiment(
            abundance_file="tests/integration/data/abundance.csv",
            glycan_file="tests/integration/data/structures.csv",
            group_file="tests/integration/data/groups.csv",
            mode="structure",
        )

    def test_run_workflow(self, experiment):
        experiment.run_workflow()

        assert isinstance(experiment.processed_abundance_table, pd.DataFrame)
        assert isinstance(experiment.meta_property_table, pd.DataFrame)
        assert isinstance(experiment.derived_trait_table, pd.DataFrame)
        assert isinstance(experiment.filtered_derived_trait_table, pd.DataFrame)
        assert isinstance(experiment.diff_results, dict)

    def test_step_by_step(self, experiment):
        # 1. Preprocess
        experiment.preprocess()
        assert isinstance(experiment.processed_abundance_table, pd.DataFrame)
        assert isinstance(experiment.meta_property_table, pd.DataFrame)

        # 2. Derive traits
        experiment.derive_traits()
        assert isinstance(experiment.derived_trait_table, pd.DataFrame)

        # 3. Filter traits
        experiment.post_filter()
        assert isinstance(experiment.filtered_derived_trait_table, pd.DataFrame)

        # 4. Differential analysis
        experiment.diff_analysis()
        assert isinstance(experiment.diff_results, dict)
