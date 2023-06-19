from collections.abc import Mapping, Callable
from pathlib import Path
from typing import Any, NoReturn, ClassVar

from .exception import ConfigError

default_config = {
    "input_file": None,
    "output_file": None,
    "filter_glycan_max_na": 0.5,
    "impute_method": "min",
    "sia_linkage": False,
    "formula_file": None,
    "filter_invalid_traits": True,
    "group_file": None,
    "structure_file": None
}

must_have = {"input_file", "output_file"}


class Config:
    """Configuration for a GlyTrait workflow.

    Valid config keys:
        - input_file: Path to the input file. Must provide.
        - output_file: Path to the output file. Must provide.
        - filter_glycan_max_na: Maximum NA percentage for a glycan to be kept.
        - impute_method: Method for imputing missing values.
        - sia_linkage: Whether to include sialic acid linkage in the analysis.
        - formula_file: Path to the formula file.
        - filter_invalid_traits: Whether to filter out invalid traits.
        - group_file: Path to the group file.
        - structure_file: Path to the structure file.
    """

    validators: ClassVar[dict[str, Callable]] = {}

    def __init__(self, user_config: Mapping[str, Any]):
        for key in must_have:
            if key not in user_config:
                raise ConfigError(f"Missing config key: {key}")
        self._config = default_config.copy()
        self.update(user_config)

    @classmethod
    def register_validator(cls, validator: Callable) -> Callable:
        """Register a validator for a config key."""
        name = validator.__name__.removeprefix("valid_")
        cls.validators[name] = validator
        return validator

    def asdict(self) -> dict[str, Any]:
        """Return the config as a dictionary."""
        return self._config.copy()

    def get(self, key: str) -> Any:
        """Get a config value."""
        return self._config[key]

    def update(self, new: Mapping[str, Any]) -> None:
        """Update the config with new values."""
        for key, value in new.items():
            if key not in self._config:
                raise KeyError(f"Invalid config key: {key}")
            self.check_validity(key, value)
        for key, value in new.items():
            self._config[key] = value

    def check_validity(self, key: str, value: Any) -> NoReturn:
        """Check if a value is valid for a config key."""
        validator = self.validators.get(key)
        if validator is None:
            return
        validator(value)

    def __repr__(self):
        return f"Config({self._config!r})"


def valid_file(
    file: str, name: str, suffix: str, *, check_exist: bool = True, none_ok: bool = True
) -> NoReturn:
    """Validate a file.

    Args:
        file: File path.
        name: Name of the file.
        suffix: File suffix.
        check_exist: Whether the file should exist.
        none_ok: Whether None is a valid value.

    Raises:
        ParameterError: If the file is invalid.
    """
    if file is None and none_ok:
        return
    if not isinstance(file, str):
        raise ConfigError(f"{name} must be a string.")
    if check_exist and not Path(file).exists():
        raise ConfigError(f"{name} does not exist.")
    if check_exist and not Path(file).is_file():
        raise ConfigError(f"{name} must be a file, not directory.")
    if Path(file).suffix != suffix:
        raise ConfigError(f"{name} must be a {suffix.upper()[1:]} file.")


@Config.register_validator
def valid_input_file(value: str) -> NoReturn:
    """Check if a value is a valid input file."""
    valid_file(value, "Input file", ".csv", check_exist=True)


@Config.register_validator
def valid_output_file(value: str) -> NoReturn:
    """Check if a value is a valid output file."""
    valid_file(value, "Output file", ".xlsx", check_exist=False)


@Config.register_validator
def valid_filter_glycan_max_na(value: float) -> NoReturn:
    """Check if a value is a valid filter_glycan_max_na."""
    if not isinstance(value, (float, int)):
        raise ConfigError("filter_glycan_max_na must be a float.")
    if not 0 <= value <= 1:
        raise ConfigError("filter_glycan_max_na must be between 0 and 1.")


@Config.register_validator
def valid_impute_method(value: str) -> NoReturn:
    """Check if a value is a valid impute_method."""
    if not isinstance(value, str):
        raise ConfigError("impute_method must be a string.")
    if value not in {"min", "mean", "median", "zero", "lod"}:
        raise ConfigError("impute_method must be one of: min, mean, median, zero, lod.")


@Config.register_validator
def valid_sia_linkage(value: bool) -> NoReturn:
    """Check if a value is a valid sia_linkage."""
    if not isinstance(value, bool):
        raise ConfigError("sia_linkage must be a boolean.")


@Config.register_validator
def valid_formula_file(value: str) -> NoReturn:
    """Check if a value is a valid formula_file."""
    valid_file(value, "Formula file", ".txt", check_exist=True)


@Config.register_validator
def valid_filter_invalid_traits(value: bool) -> NoReturn:
    """Check if a value is a valid filter_invalid_traits."""
    if not isinstance(value, bool):
        raise ConfigError("filter_invalid_traits must be a boolean.")


@Config.register_validator
def valid_group_file(value: str) -> NoReturn:
    """Check if a value is a valid group_file."""
    valid_file(value, "Group file", ".csv", check_exist=True)


@Config.register_validator
def valid_structure_file(value: str) -> NoReturn:
    """Check if a value is a valid structure_file."""
    valid_file(value, "Structure file", ".csv", check_exist=True)
