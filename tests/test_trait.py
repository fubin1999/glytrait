import numpy as np
import pandas as pd
import pytest

from glytrait import trait
from glytrait.exception import *
from . import glycoct


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
        with pytest.raises(FormulaError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=["invalid", "isComplex"],
                denominator_properties=["isComplex"],
            )
        msg = "`numerator_properties` contains invalid meta properties: invalid."
        assert msg in str(excinfo.value)

    def test_init_0_length_properties(self):
        with pytest.raises(FormulaError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=[],
                denominator_properties=["isHybrid"],
            )
        assert "`numerator_properties` cannot be empty." in str(excinfo.value)

    def test_init_dot_in_numerator(self):
        with pytest.raises(FormulaError) as excinfo:
            trait.TraitFormula(
                description="The ratio of high-mannose to hybrid glycans",
                name="MHy",
                numerator_properties=["."],
                denominator_properties=["isHybrid"],
            )
        assert "'.' should not be used in the numerator." in str(excinfo.value)

    def test_init_dot_with_others_in_denominator(self):
        with pytest.raises(FormulaError) as excinfo:
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

    def test_calcu_trait_inf(self, formula1, meta_property_table, abundance_table):
        meta_property_table = meta_property_table.drop("G4", axis=0)
        abundance_table = abundance_table.drop("G4", axis=1)
        formula1.initialize(meta_property_table)
        result = formula1.calcu_trait(abundance_table)
        expected = np.array([np.nan, np.nan, np.nan])
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
    with pytest.raises(FormulaError) as excinfo:
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
    with pytest.raises(FormulaError):
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


def test_load_formulas_with_user_file(clean_dir):
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


def test_load_formulas_without_user_file():
    result = list(trait.load_formulas(None))
    assert len(result) == 164


def test_load_formulas_bad_formula(clean_dir):
    content = """@ Relative abundance of complex type glycans within total spectrum
$ TC = (isComplex) / (. * isComplex)
"""
    user_file = clean_dir / "formulas.txt"
    user_file.write_text(content)
    with pytest.raises(FormulaError) as excinfo:
        list(trait.load_formulas(user_file))
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
        list(trait.load_formulas(user_file))
    assert "Invalid line: '@ A duplicate of TM'" in str(excinfo.value)
    assert "One description line must follow a formula line." in str(excinfo.value)


def test_save_trait_formula_template(clean_dir):
    trait.save_trait_formula_template(clean_dir)
    template_file = clean_dir / "trait_formula.txt"
    assert template_file.exists()
    assert "Trait Formula Overview" in template_file.read_text()


def test_calcu_trait():
    trait_formulas = [
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

    abundance_df = pd.DataFrame(
        {
            "G1": [1, 2, 3],
            "G2": [4, 5, 6],
            "G3": [7, 8, 9],
        },
        index=["S1", "S2", "S3"],
    )

    meta_prop_df = pd.DataFrame(
        {
            "isHighMannose": [True, False, False],
            "isHybrid": [False, False, True],
            "isComplex": [False, True, False],
        },
        index=["G1", "G2", "G3"],
    )

    result = trait.calcu_derived_trait(abundance_df, meta_prop_df, trait_formulas)
    expected = pd.DataFrame(
        {
            "TM": [1 / 12, 2 / 15, 3 / 18],
            "MHy": [1 / 7, 2 / 8, 3 / 9],
        },
        index=["S1", "S2", "S3"],
    )
    pd.testing.assert_frame_equal(result, expected)


def test_build_meta_property_table_no_sia_linkage(make_glycan):
    glycoct_strings = {
        k: v for k, v in glycoct.__dict__.items() if k.startswith("test_glycoct")
    }
    del glycoct_strings["test_glycoct_13"]  # an O-glycan

    glycan_ids = list(glycoct_strings.keys())
    glycans = [make_glycan(s) for s in glycoct_strings.values()]

    result = trait.build_meta_property_table(glycan_ids, glycans, sia_linkage=False)
    expected = pd.DataFrame(
        {
            "isComplex": [1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1],
            "isHighMannose": [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            "isHybrid": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            "isBisecting": [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            "is1Antennary": [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            "is2Antennary": [1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            "is3Antennary": [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            "is4Antennary": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            "totalAntenna": [2, 2, 0, 0, 0, 1, 0, 3, 2, 2, 2, 2],
            "coreFuc": [0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
            "antennaryFuc": [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
            "hasAntennaryFuc": [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            "totalFuc": [0, 1, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0],
            "hasFuc": [0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
            "noFuc": [1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1],
            "totalSia": [2, 1, 0, 1, 0, 1, 0, 1, 1, 2, 2, 2],
            "hasSia": [1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1],
            "noSia": [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
            "totalMan": [3, 3, 5, 5, 3, 3, 6, 3, 3, 3, 3, 3],
            "totalGal": [2, 2, 0, 1, 0, 1, 0, 3, 1, 2, 2, 2],
        },
        index=glycan_ids,
    )
    pd.testing.assert_frame_equal(result, expected, check_dtype=False)


def test_build_meta_property_table_sia_linkage(make_glycan):
    glycoct_strings = {
        k: v
        for k, v in glycoct.__dict__.items()
        if k in ("test_glycoct_10", "test_glycoct_11", "test_glycoct_12")
    }

    glycan_ids = list(glycoct_strings.keys())
    glycans = [make_glycan(s) for s in glycoct_strings.values()]

    result = trait.build_meta_property_table(glycan_ids, glycans, sia_linkage=True)
    partial_expected = pd.DataFrame(
        {
            "a23Sia": [0, 2, 1],
            "a26Sia": [2, 0, 1],
            "hasa23Sia": [0, 1, 1],
            "hasa26Sia": [1, 0, 1],
            "noa23Sia": [1, 0, 0],
            "noa26Sia": [0, 1, 0],
        },
        index=glycan_ids,
    )
    assert len(result.columns) == 26
    partial_result = result[partial_expected.columns]
    pd.testing.assert_frame_equal(partial_result, partial_expected, check_dtype=False)


def test_calcu_direct_trait():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3],
            "B": [4, 5, 6],
            "C": [7, 8, 9],
        }
    )
    result = trait.calcu_direct_trait(df)
    expected = pd.DataFrame(
        {
            "A": [1 / 12, 2 / 15, 3 / 18],
            "B": [4 / 12, 5 / 15, 6 / 18],
            "C": [7 / 12, 8 / 15, 9 / 18],
        }
    )
    pd.testing.assert_frame_equal(result, expected, check_dtype=False)


def test_filter_derived_trait():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3],
            "B": [np.nan, np.nan, np.nan],
            "C": [1, 1, np.nan],
            "D": [0, 0, 0],
            "E": [1, 1, 1],
            "F": [0.5, 0.5, 0.5],
        }
    )
    result = trait.filter_derived_trait(df)
    expected = pd.DataFrame(
        {
            "A": [1, 2, 3],
        }
    )
    pd.testing.assert_frame_equal(result, expected, check_dtype=False)
