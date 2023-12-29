import pytest

from glytrait.api import _Config


class TestConfig:
    @pytest.fixture
    def config_dict(self) -> dict:
        return dict(
            mode="structure",
            filter_max_na=0.0,
            impute_method="zero",
            post_filtering=True,
            correlation_threshold=0.9,
            sia_linkage=False,
            custom_formula_file=None,
        )

    def test_init(self, config_dict):
        config = _Config(**config_dict)
        assert config.mode == "structure"
        assert config.filter_max_na == 0.0
        assert config.impute_method == "zero"
        assert config.post_filtering is True
        assert config.correlation_threshold == 0.9
        assert config.sia_linkage is False
        assert config.custom_formula_file is None

    @pytest.mark.parametrize("mode", ["invalid", 1, None])
    def test_invalid_mode(self, config_dict, mode):
        config_dict["mode"] = mode
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("filter_max_na", ["invalid", -1, 2, None])
    def test_invalid_filter_max_na(self, config_dict, filter_max_na):
        config_dict["filter_max_na"] = filter_max_na
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("impute_method", ["invalid", 1, None])
    def test_invalid_impute_method(self, config_dict, impute_method):
        config_dict["impute_method"] = impute_method
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("post_filtering", ["invalid", 1, None])
    def test_invalid_post_filtering(self, config_dict, post_filtering):
        config_dict["post_filtering"] = post_filtering
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("correlation_threshold", ["invalid", -2, -0.5, 2, None])
    def test_invalid_correlation_threshold(self, config_dict, correlation_threshold):
        config_dict["correlation_threshold"] = correlation_threshold
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("sia_linkage", ["invalid", 1, None])
    def test_invalid_sia_linkage(self, config_dict, sia_linkage):
        config_dict["sia_linkage"] = sia_linkage
        with pytest.raises(ValueError):
            _Config(**config_dict)

    @pytest.mark.parametrize("custom_formula_file", [1])
    def test_invalid_custom_formula_file(self, config_dict, custom_formula_file):
        config_dict["custom_formula_file"] = custom_formula_file
        with pytest.raises(ValueError):
            _Config(**config_dict)
