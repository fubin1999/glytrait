import numpy as np
import pandas as pd
import pytest
from attrs import evolve

import glytrait.formula as fml
from glytrait.exception import FormulaError


@pytest.fixture
def write_content(clean_dir):
    def _write_content(content):
        file = clean_dir / "formulas.txt"
        file.write_text(content)
        return file

    return _write_content


@pytest.mark.skip("`TraitFormula` to be updated.")
class TestTraitFormula:
    @pytest.fixture
    def formula1(self):
        return fml.TraitFormula(
            description="The ratio of high-mannose to hybrid glycans",
            name="MHy",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isHybrid"],
        )

    @pytest.fixture
    def formula2(self):
        return fml.TraitFormula(
            description="Relative abundance of high mannose type glycans within total spectrum",
            name="TM",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["."],
        )

    @pytest.fixture
    def meta_property_table(self):
        data = {
            ".": [1, 1, 1, 1],
            "isComlex": [False, False, True, False],
            "isHighMannose": [True, True, False, False],
            "isHybrid": [False, False, False, True],
        }
        return pd.DataFrame(data, index=["G1", "G2", "G3", "G4"])

    @pytest.fixture
    def abundance_table(self):
        data = dict(
            G1=[1, 2, 3],
            G2=[4, 5, 6],
            G3=[7, 8, 9],
            G4=[10, 11, 12],
        )
        return pd.DataFrame(data, index=["S1", "S2", "S3"], dtype=float)

    def test_init_invalid_properties(self):
        with pytest.raises(FormulaError) as excinfo:
            fml.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=["invalid", "isComplex"],
                denominator_properties=["isComplex"],
            )
        msg = "`numerator_properties` contains invalid meta properties: invalid."
        assert msg in str(excinfo.value)

    def test_init_0_length_properties(self):
        with pytest.raises(FormulaError) as excinfo:
            fml.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=[],
                denominator_properties=["isHybrid"],
            )
        assert "`numerator_properties` cannot be empty." in str(excinfo.value)

    def test_init_dot_in_numerator(self):
        with pytest.raises(FormulaError) as excinfo:
            fml.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=["."],
                denominator_properties=["isHybrid"],
            )
        assert "'.' should not be used in the numerator." in str(excinfo.value)

    def test_init_dot_with_others_in_denominator(self):
        with pytest.raises(FormulaError) as excinfo:
            fml.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
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
            fml.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=["isHighMannose"],
                denominator_properties=["isHybrid"],
                coefficient=coef,
            )

    def test_init_invalid_type(self):
        with pytest.raises(ValueError):
            fml.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="invalid",
                numerator_properties=["isHighMannose"],
                denominator_properties=["isHybrid"],
            )

    def test_init_wrong_type(self):
        with pytest.raises(FormulaError) as excinfo:
            fml.TraitFormula(
                description="Should be a composition trait",
                name="SomeTrait",
                type="composition",
                numerator_properties=["isHighMannose"],
                denominator_properties=["."],
            )
        msg = "`numerator_properties` contains invalid meta properties: isHighMannose."
        assert msg in str(excinfo.value)

    @pytest.mark.parametrize(
        "type, numerator, denominator, expected",
        [
            ("structure", ["isHighMannose"], ["isHybrid"], False),
            ("structure", ["hasa23Sia"], ["."], True),
            ("composition", ["isHighBranching"], ["."], False),
            ("composition", ["hasa23Sia"], ["."], True),
        ],
    )
    def test_sia_linkage(self, type, numerator, denominator, expected):
        formula = fml.TraitFormula(
            description="Should be a composition trait",
            name="SomeTrait",
            type=type,
            numerator_properties=numerator,
            denominator_properties=denominator,
        )
        assert formula.sia_linkage == expected

    def test_calcu_trait_without_initialization(self, formula1):
        with pytest.raises(RuntimeError):
            formula1.calcu_trait(None)

    def test_calcu_trait(self, formula1, meta_property_table, abundance_table):
        formula1.initialize(meta_property_table)
        result = formula1.calcu_trait(abundance_table)
        expected = [5 / 10, 7 / 11, 9 / 12]
        np.testing.assert_array_equal(result, expected)

    def test_calcu_trait_glycan_order_not_same(
        self, formula1, meta_property_table, abundance_table
    ):
        meta_property_table = meta_property_table.reindex(["G2", "G1", "G3", "G4"])
        formula1.initialize(meta_property_table)
        with pytest.raises(AssertionError) as excinfo:
            formula1.calcu_trait(abundance_table)

    def test_calcu_trait_with_dot(self, formula2, meta_property_table, abundance_table):
        formula2.initialize(meta_property_table)
        result = formula2.calcu_trait(abundance_table)
        expected = [5 / 22, 7 / 26, 9 / 30]
        np.testing.assert_array_equal(result, expected)

    def test_calcu_trait_with_coef(
        self, formula1, meta_property_table, abundance_table
    ):
        formula = evolve(formula1, coefficient=2)
        formula.initialize(meta_property_table)
        result = formula.calcu_trait(abundance_table)
        expected = [10 / 10, 14 / 11, 18 / 12]
        np.testing.assert_array_equal(result, expected)

    def test_calcu_trait_inf(self, formula1, meta_property_table, abundance_table):
        meta_property_table = meta_property_table.drop("G4", axis=0)
        abundance_table = abundance_table.drop("G4", axis=1)
        formula1.initialize(meta_property_table)
        result = formula1.calcu_trait(abundance_table)
        expected = np.array([np.nan, np.nan, np.nan])
        np.testing.assert_array_equal(result, expected)

    def test_try_to_change_numerator(self, formula1):
        numerator = formula1.numerator_properties
        before_change = numerator.copy()
        numerator.append("new_property")
        assert formula1.numerator_properties == before_change


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
    (
        r_name,
        r_num_props,
        r_den_props,
        r_coef,
    ) = fml.parse_formula_expression(expression)
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
    with pytest.raises(FormulaError) as excinfo:
        fml.parse_formula_expression(expression)
    assert f"Invalid expression: '{expression}'" in str(excinfo.value)


@pytest.mark.skip("`TraitFormula` to be updated.")
def test_load_default_formulas():
    structure_formulas = list(fml.load_default_formulas("structure"))
    composition_formulas = list(fml.load_default_formulas("composition"))
    assert len(structure_formulas) > 0
    assert len(composition_formulas) > 0
    assert structure_formulas[0].type == "structure"
    assert composition_formulas[0].type == "composition"
    assert len(structure_formulas) != len(composition_formulas)


@pytest.mark.skip("`TraitFormula` to be updated.")
class TestLoadFormulasFromFile:
    """Test `load_formulas_from_file` function."""

    def test_basic(self, write_content):
        description = "The ratio of high-mannose to hybrid glycans"
        expression = "MHy = (isHighMannose) / (isHybrid)"
        file = write_content(f"@ {description}\n$ {expression}")
        result = list(fml.load_formulas_from_file(file, "structure"))
        assert len(result) == 1
        assert result[0].description == description
        assert result[0].name == "MHy"

    def test_duplicated_formulas(self, write_content):
        description1 = "The ratio of high-mannose to hybrid glycans"
        expression1 = "MHy = (isHighMannose) / (isHybrid)"
        description2 = (
            "Relative abundance of high mannose type glycans within total spectrum"
        )
        expression2 = "MHy = (isHighMannose) / (isHybrid)"  # same expression
        content = (
            f"@ {description1}\n$ {expression1}\n@ {description2}\n$ {expression2}"
        )
        file = write_content(content)
        with pytest.raises(FormulaError) as excinfo:
            list(fml.load_formulas_from_file(file, "structure"))
        assert "Duplicate formula name: MHy." in str(excinfo.value)


def test_save_builtin_formula(clean_dir):
    fml.save_builtin_formula(clean_dir)
    struc_file = clean_dir / "struc_builtin_formulas.txt"
    comp_file = clean_dir / "comp_builtin_formulas.txt"
    assert struc_file.exists()
    assert comp_file.exists()


class TestDeconvoluteFormulaFile:
    """Test `deconvolute_formula_file` function."""

    def test_basic(self, write_content):
        file = write_content("@ Description\n$ Expression\n")
        result = list(fml.deconvolute_formula_file(file))
        expected = [("Description", "Expression")]
        assert result == expected

    def test_first_line_not_description(self, write_content):
        file = write_content("$ Expression\n@ Description\n")
        with pytest.raises(FormulaError) as excinfo:
            list(fml.deconvolute_formula_file(file))
        assert "No description before expression 'Expression'" in str(excinfo.value)

    def test_two_descriptions(self, write_content):
        file = write_content("@ Description1\n@ Description2\n")
        with pytest.raises(FormulaError) as excinfo:
            list(fml.deconvolute_formula_file(file))
        assert "No expression follows description 'Description1'." in str(excinfo.value)

    def test_two_expressions(self, write_content):
        file = write_content("@ Description\n$ Expression1\n$ Expression2")
        with pytest.raises(FormulaError) as excinfo:
            list(fml.deconvolute_formula_file(file))
        assert "No description before expression 'Expression2'." in str(excinfo.value)

    def test_no_last_expression(self, write_content):
        file = write_content("@ Description1\n$ Expression1\n@Description2")
        with pytest.raises(FormulaError) as excinfo:
            list(fml.deconvolute_formula_file(file))
        assert "No expression follows description 'Description2'." in str(excinfo.value)
