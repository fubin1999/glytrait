import numpy as np
import pandas as pd
import pytest
from attrs import evolve

import glytrait.formula
from glytrait.exception import FormulaError


class TestTraitFormula:
    @pytest.fixture
    def formula1(self):
        return glytrait.formula.TraitFormula(
            description="The ratio of high-mannose to hybrid glycans",
            name="MHy",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isHybrid"],
        )

    @pytest.fixture
    def formula2(self):
        return glytrait.formula.TraitFormula(
            description="Relative abundance of high mannose type glycans within total spectrum",
            name="TM",
            type="structure",
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
        with pytest.raises(FormulaError) as excinfo:
            glytrait.formula.TraitFormula(
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
            glytrait.formula.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=[],
                denominator_properties=["isHybrid"],
            )
        assert "`numerator_properties` cannot be empty." in str(excinfo.value)

    def test_init_dot_in_numerator(self):
        with pytest.raises(FormulaError) as excinfo:
            glytrait.formula.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=["."],
                denominator_properties=["isHybrid"],
            )
        assert "'.' should not be used in the numerator." in str(excinfo.value)

    def test_init_dot_with_others_in_denominator(self):
        with pytest.raises(FormulaError) as excinfo:
            glytrait.formula.TraitFormula(
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
            glytrait.formula.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="structure",
                numerator_properties=["isHighMannose"],
                denominator_properties=["isHybrid"],
                coefficient=coef,
            )

    def test_init_invalid_type(self):
        with pytest.raises(ValueError):
            glytrait.formula.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                type="invalid",
                numerator_properties=["isHighMannose"],
                denominator_properties=["isHybrid"],
            )

    def test_init_wrong_type(self):
        with pytest.raises(FormulaError) as excinfo:
            glytrait.formula.TraitFormula(
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
        formula = glytrait.formula.TraitFormula(
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
        "trait1, trait2, expected",
        [
            ("A2G", "CG", True),
            ("A2F", "CF", True),
            ("A2Fa", "CFa", True),
            ("A2Fc", "CFc", True),
            ("CFa", "CF", True),
            ("A2S", "CS", True),
            ("A2E", "A2S", True),
            ("A2E", "CE", True),
            ("A2SG", "A2G", True),
            ("A2FSG", "A2FG", True),
            ("A2FSG", "A2SG", True),
            ("A2F0G", "A2FG", False),
            ("A2FSG", "A2G", False),
        ],
    )
    def test_is_child_of(self, trait1, trait2, expected):
        formulas = glytrait.formula.load_formulas("structure")
        formula_map = {f.name: f for f in formulas}
        formulas1 = formula_map[trait1]
        formulas2 = formula_map[trait2]
        assert formulas1.is_child_of(formulas2) == expected


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
    r_name, r_num_props, r_den_props, r_coef = glytrait.formula._parse_expression(
        expression
    )
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
        glytrait.formula._parse_expression(expression)
    assert f"Invalid expression: '{expression}'" in str(excinfo.value)


@pytest.mark.parametrize(
    "type, expected",
    [
        ("structure", 299),
        ("composition", 75),
    ],
)
def test_load_default_formulas(type, expected):
    result = list(glytrait.formula._load_default_formulas(type=type))
    assert len(result) == expected


def test_load_user_formulas(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (.)

@ A duplicate of the above
$ TC = (isComplex) / (isHybrid)"""
    file = clean_dir / "formulas.txt"
    file.write_text(content)
    result = list(glytrait.formula._load_user_formulas(file, type="structure"))
    assert len(result) == 1
    assert result[0].name == "TC"
    assert (
        result[0].description
        == "Relative abundance of complex type glycans within total spectrum"
    )


def test_load_formulas_with_user_file(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (.)

@ A duplicate of TM
$ TM = (isHighMannose) / (.)"""
    user_file = clean_dir / "formulas.txt"
    user_file.write_text(content)
    result = list(glytrait.formula.load_formulas("structure", user_file))

    TM = [f for f in result if f.name == "TM"]
    assert len(TM) == 1
    assert (
        TM[0].description
        == "Relative abundance of high mannose type glycans within total spectrum"
    )

    assert "TC" in [f.name for f in result]


def test_load_formulas_without_user_file():
    result = list(glytrait.formula.load_formulas(type="structure"))
    assert len(result) == 299


def test_load_formulas_bad_formula(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (. * isComplex)
"""
    user_file = clean_dir / "formulas.txt"
    user_file.write_text(content)
    with pytest.raises(FormulaError) as excinfo:
        list(glytrait.formula.load_formulas("structure", user_file))
    assert "Invalid line: '$ TC = (isComplex) / (. * isComplex)'" in str(excinfo.value)
    assert "Error in formula 'TC = (isComplex) / (. * isComplex)'" in str(excinfo.value)


def test_load_formulas_end_with_description(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (.)
    
@ A duplicate of TM
"""
    user_file = clean_dir / "formulas.txt"
    user_file.write_text(content)
    with pytest.raises(FormulaError) as excinfo:
        list(glytrait.formula.load_formulas("structure", user_file))
    assert "Invalid line: '@ A duplicate of TM'" in str(excinfo.value)
    assert "One description line must follow a formula line." in str(excinfo.value)


def test_save_trait_formula_template(clean_dir):
    glytrait.formula.save_trait_formula_template(clean_dir)
    template_file = clean_dir / "trait_formula.txt"
    assert template_file.exists()
    assert "Trait Formula Overview" in template_file.read_text()


def test_save_builtin_formula(clean_dir):
    glytrait.formula.save_builtin_formula(clean_dir)
    struc_file = clean_dir / "struc_builtin_formulas.txt"
    comp_file = clean_dir / "comp_builtin_formulas.txt"
    assert struc_file.exists()
    assert comp_file.exists()
