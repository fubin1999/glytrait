import numpy as np
import pandas as pd
import pytest

from glytrait import trait


class TestTraitFormula:
    @pytest.fixture
    def formula1(self):
        return trait.TraitFormula(
            description="The ratio of high-mannose to hybrid glycans",
            name="MHy",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isHybrid"],
        )

    @pytest.fixture
    def formula2(self):
        return trait.TraitFormula(
            description="Relative abundance of high mannose type glycans within total spectrum",
            name="TM",
            numerator_properties=["isHighMannose"],
            denominator_properties=["."],
        )

    @pytest.fixture
    def meta_property_table(self):
        data = dict(
            isComlex=[False, False, True, False],
            isHighMannose=[True, True, False, False],
            isHybrid=[False, False, False, True],
        )
        return pd.DataFrame(data, index=["G1", "G2", "G3", "G4"])

    @pytest.fixture
    def abundance_table(self):
        data = dict(
            G1=[1, 2, 3],
            G2=[4, 5, 6],
            G3=[7, 8, 9],
            G4=[10, 11, 12],
        )
        return pd.DataFrame(data, index=["S1", "S2", "S3"])

    def test_init_invalid_properties(self):
        with pytest.raises(ValueError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=["invalid", "isComplex"],
                denominator_properties=["isComplex"],
            )
        msg = "`numerator_properties` contains invalid meta properties: invalid."
        assert msg in str(excinfo.value)

    def test_init_0_length_properties(self):
        with pytest.raises(ValueError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=[],
                denominator_properties=["isHybrid"],
            )
        assert "`numerator_properties` cannot be empty." in str(excinfo.value)

    def test_calcu_trait_without_initialization(self, formula1):
        with pytest.raises(RuntimeError):
            formula1.calcu_trait(None)

    def test_calcu_trait(self, formula1, meta_property_table, abundance_table):
        formula1.initialize(meta_property_table)
        result = formula1.calcu_trait(abundance_table)
        expected = [5 / 10, 7 / 11, 9 / 12]
        np.testing.assert_array_equal(result, expected)

    def test_calcu_trait_with_dot(self, formula2, meta_property_table, abundance_table):
        formula2.initialize(meta_property_table)
        result = formula2.calcu_trait(abundance_table)
        expected = [5 / 22, 7 / 26, 9 / 30]
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "expression, name, num_props, den_props",
    [
        ("TM = (isHighMannose) / (.)", "TM", ["isHighMannose"], ["."]),
        ("MHy = (isHighMannose) / (isHybrid)", "MHy", ["isHighMannose"], ["isHybrid"]),
        (
            "CA1 = (isComplex * is1Antennay) / (isComplex)",
            "CA1",
            ["isComplex", "is1Antennay"],
            ["isComplex"],
        ),
    ],
)
def test_parse_expression(expression, name, num_props, den_props):
    result = trait._parse_expression(expression)
    expected = (name, num_props, den_props)
    assert result == expected


@pytest.mark.parametrize(
    "expression",
    [
        "(isHighMannose) / (.)",
        "TM = (isHighMannose) / (isHybrid",
        "MHy = (isHighMannose)  (isHybrid)",
        "CA1 = (isComplex  is1Antennay) / (isComplex)",
        "CA1 * TM = (isComplex * is1Antennay) / (isComplex)",
    ],
)
def test_parse_expression_invalid(expression):
    with pytest.raises(trait.FormulaParseError) as excinfo:
        trait._parse_expression(expression)
    assert f"Invalid expression: '{expression}'" in str(excinfo.value)


def test_load_default_formulas(tmp_path, monkeypatch):
    text = """
# Some comment

@ Relative abundance of high mannose type glycans within total spectrum
$ TM = (isHighMannose) / (.)

@ The ratio of high-mannose to hybrid glycans
$ MHy = (isHighMannose) / (isHybrid)
"""
    formula_file = tmp_path / "formula.txt"
    formula_file.write_text(text)
    monkeypatch.setattr(trait, "DEFAULT_FORMULA_FILE", str(formula_file))
    result = list(trait.load_default_formulas())
    expected = [
        trait.TraitFormula(
            description="Relative abundance of high mannose type glycans within total spectrum",
            name="TM",
            numerator_properties=["isHighMannose"],
            denominator_properties=["."],
        ),
        trait.TraitFormula(
            description="The ratio of high-mannose to hybrid glycans",
            name="MHy",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isHybrid"],
        ),
    ]
    assert result == expected


@pytest.mark.parametrize(
    "text",
    [
        """
@ Description 1
@ Description 2
""",
        """
@ Description 1
$ TM = (isHighMannose) / (.)
$ MHy = (isHighMannose) / (isHybrid)
""",
    ],
)
def test_load_default_formulas_invalid(text, tmp_path, monkeypatch):
    formula_file = tmp_path / "formula.txt"
    formula_file.write_text(text)
    monkeypatch.setattr(trait, "DEFAULT_FORMULA_FILE", str(formula_file))
    with pytest.raises(trait.FormulaParseError):
        list(trait.load_default_formulas())
