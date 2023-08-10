import pytest

import glytrait.workflow
from glytrait.config import default_config as glytrait_default_config, Config


@pytest.fixture
def default_config(monkeypatch) -> Config:
    monkeypatch.setattr(glytrait.workflow.Config, "validators", [])
    new = dict(
        input_file="input_file",
        output_file="output_file",
        mode="structure",
        filter_glycan_max_na=0.5,
        impute_method="min",
        corr_threshold=1.0,
        corr_method="pearson",
        sia_linkage=False,
        formula_file=None,
        post_filtering=True,
        group_file=None,
        structure_file=None,
        database=None,
    )
    config_dict = glytrait_default_config.copy()
    config_dict.update(new)
    config = Config(config_dict)
    return config
