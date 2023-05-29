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

    def test_init_dot_in_numerator(self):
        with pytest.raises(ValueError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=["."],
                denominator_properties=["isHybrid"],
            )
        assert "'.' should not be used in the numerator." in str(excinfo.value)

    def test_init_dot_with_others_in_denominator(self):
        with pytest.raises(ValueError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=["isHighMannose"],
                denominator_properties=["isHybrid", "."],
            )
        assert (
            "'.' should not be used with other meta properties in the denominator."
            in str(excinfo.value)
        )

    @pytest.mark.parametrize("coef", [-1, 0])
    def test_init_invalid_coef(self, coef):
        with pytest.raises(ValueError):
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=["isHighMannose"],
                denominator_properties=["isHybrid"],
                coefficient=coef,
            )

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

    def test_calcu_trait_with_coef(
        self, formula1, meta_property_table, abundance_table
    ):
        formula1.coefficient = 2
        formula1.initialize(meta_property_table)
        result = formula1.calcu_trait(abundance_table)
        expected = [10 / 10, 14 / 11, 18 / 12]
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "expression, name, num_props, den_props, coef",
    [
        ("TM = (isHighMannose) / (.)", "TM", ["isHighMannose"], ["."], 1.0),
        (
            "MHy = (isHighMannose) / (isHybrid)",
            "MHy",
            ["isHighMannose"],
            ["isHybrid"],
            1.0,
        ),
        (
            "CA1 = (isComplex * is1Antennay) / (isComplex)",
            "CA1",
            ["isComplex", "is1Antennay"],
            ["isComplex"],
            1.0,
        ),
        (
            "CA1 = (is1Antennay) // (isComplex)",
            "CA1",
            ["isComplex", "is1Antennay"],
            ["isComplex"],
            1.0,
        ),
        (
            "CA1 = (isComplex * is1Antennay) / (isComplex) * 1/2",
            "CA1",
            ["isComplex", "is1Antennay"],
            ["isComplex"],
            0.5,
        ),
        (
            "CA1 = (isComplex * is1Antennay) / (isComplex) * 0.5",
            "CA1",
            ["isComplex", "is1Antennay"],
            ["isComplex"],
            0.5,
        ),
        (
            "CA1 = (isComplex * is1Antennay) / (isComplex) * 2",
            "CA1",
            ["isComplex", "is1Antennay"],
            ["isComplex"],
            2,
        ),
    ],
)
def test_parse_expression(expression, name, num_props, den_props, coef):
    r_name, r_num_props, r_den_props, r_coef = trait._parse_expression(expression)
    assert r_name == name
    assert sorted(r_num_props) == sorted(num_props)
    assert sorted(r_den_props) == sorted(den_props)
    assert pytest.approx(coef) == r_coef


@pytest.mark.parametrize(
    "expression",
    [
        "(isHighMannose) / (.)",
        "TM = (isHighMannose) / (isHybrid",
        "MHy = (isHighMannose)  (isHybrid)",
        "CA1 = (isComplex  is1Antennay) / (isComplex)",
        "CA1 * TM = (isComplex * is1Antennay) / (isComplex)",
        "MHy = (isHighMannose) / (isHybrid) * a",
        "MHy = (isHighMannose) / (isHybrid) * ",
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


def test_load_user_formulas(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (.)

@ A duplicate of the above
$ TC = (isComplex) / (isHybrid)"""
    file = clean_dir / "formulas.txt"
    file.write_text(content)
    result = list(trait.load_user_formulas(file))
    assert len(result) == 1
    assert result[0].name == "TC"
    assert (
        result[0].description
        == "Relative abundance of complex type glycans within total spectrum"
    )


def test_load_formulas(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (.)

@ A duplicate of TM
$ TM = (isHighMannose) / (.)"""
    user_file = clean_dir / "formulas.txt"
    user_file.write_text(content)
    result = list(trait.load_formulas(user_file))

    TM = [f for f in result if f.name == "TM"]
    assert len(TM) == 1
    assert (
        TM[0].description
        == "Relative abundance of high mannose type glycans within total spectrum"
    )

    assert "TC" in [f.name for f in result]
