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


@pytest.fixture
def mp_table() -> pd.DataFrame:
    data = {
        "mp_int": [1, 2, 3],
        "mp_bool": [True, False, True],
        "mp_str": ["a", "b", "c"],
    }
    df = pd.DataFrame(data, index=["G1", "G2", "G3"])
    df = df.astype({"mp_int": "UInt8", "mp_bool": "boolean", "mp_str": "category"})
    return df


class TestConstantTerm:

    def test_expr(self):
        term = fml.ConstantTerm(1)
        assert term.expr == "1"

    def test_call(self, mp_table):
        term = fml.ConstantTerm(1)
        result = term(mp_table)
        expected = pd.Series([1, 1, 1], index=mp_table.index, name="1", dtype="UInt8")
        pd.testing.assert_series_equal(result, expected)


class TestNumericalTerm:

    def test_expr(self):
        term = fml.NumericalTerm("mp_int")
        assert term.expr == "mp_int"

    def test_call(self, mp_table):
        term = fml.NumericalTerm("mp_int")
        result = term(mp_table)
        expected = pd.Series(
            [1, 2, 3], index=mp_table.index, name="mp_int", dtype="UInt8"
        )
        pd.testing.assert_series_equal(result, expected)

    @pytest.mark.parametrize("mp", ["mp_bool", "mp_str"])
    def test_call_wrong_type(self, mp_table, mp):
        term = fml.NumericalTerm(mp)
        with pytest.raises(fml.FormulaError):
            term(mp_table)

    def test_call_mp_not_exist(self, mp_table):
        term = fml.NumericalTerm("mp_not_exist")
        with pytest.raises(fml.FormulaError):
            term(mp_table)


class TestCompareTerm:

    @pytest.mark.parametrize(
        "mp, operator, value, expected",
        [
            ("mp_int", ">", 2, "mp_int > 2"),
            ("mp_bool", "==", True, "mp_bool == True"),
            ("mp_str", "!=", "b", "mp_str != 'b'"),
        ],
    )
    def test_expr(self, mp, operator, value, expected):
        term = fml.CompareTerm(mp, operator, value)
        assert term.expr == expected

    @pytest.mark.parametrize(
        "mp, operator, value, expected",
        [
            ("mp_int", ">", 2, [0, 0, 1]),
            ("mp_int", ">=", 2, [0, 1, 1]),
            ("mp_int", "<", 2, [1, 0, 0]),
            ("mp_int", "<=", 2, [1, 1, 0]),
            ("mp_int", "==", 2, [0, 1, 0]),
            ("mp_int", "!=", 2, [1, 0, 1]),
            ("mp_bool", "==", True, [1, 0, 1]),
            ("mp_bool", "!=", True, [0, 1, 0]),
            ("mp_str", "==", "b", [0, 1, 0]),
            ("mp_str", "!=", "b", [1, 0, 1]),
        ],
    )
    def test_call(self, mp_table, mp, operator, value, expected):
        term = fml.CompareTerm(mp, operator, value)
        result = term(mp_table)
        expected = pd.Series(
            expected, index=mp_table.index, name=term.expr, dtype="UInt8"
        )
        pd.testing.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "mp, operator, value",
        [
            ("mp_bool", ">", 2),
            ("mp_bool", ">=", 2),
            ("mp_bool", "<", 2),
            ("mp_bool", "<=", 2),
            ("mp_str", ">", 2),
            ("mp_str", ">=", 2),
            ("mp_str", "<", 2),
            ("mp_str", "<=", 2),
        ],
    )
    def test_type_error(self, mp_table, mp, operator, value):
        term = fml.CompareTerm(mp, operator, value)
        with pytest.raises(fml.FormulaError):
            term(mp_table)


class TestParseFormulaExpression:
    """Test `_parse_formula_expression` function."""

    @pytest.mark.parametrize(
        "expr, numerators, denominators",
        [
            (
                "A = (mp1 == 1) // (mp2 > 2)",
                ["mp1 == 1", "mp2 > 2"],
                ["mp2 > 2"],
            ),
            (
                "A = (mp1 == 1) / (mp2 > 2)",
                ["mp1 == 1"],
                ["mp2 > 2"],
            ),
            (
                "A = (mp1 == 1) / 1",
                ["mp1 == 1"],
                ["1"],
            ),
            (
                "A = (mp1 == 1) / (mp2 == 'b')",
                ["mp1 == 1"],
                ["mp2 == 'b'"],
            ),
            (
                'A = (mp1 == 1) / (mp2 == "b")',
                ["mp1 == 1"],
                ["mp2 == 'b'"],
            ),
            (
                "A = (mp1 == 1) / (mp2 == True)",
                ["mp1 == 1"],
                ["mp2 == True"],
            ),
            (
                "A = (mp1 == 1) * (mp2 > 2)/ 1",
                ["mp1 == 1", "mp2 > 2"],
                ["1"],
            ),
            (
                "A = (mp1 == 1) * (mp2 > 2) * (mp3 <= 4) / 1",
                ["mp1 == 1", "mp2 > 2", "mp3 <= 4"],
                ["1"],
            ),
            (
                "A = (mp1 == 1) * 2 / 1",
                ["mp1 == 1", "2"],
                ["1"],
            ),
            (
                "A = mp1 / mp2",
                ["mp1"],
                ["mp2"],
            ),
            (
                "A = (mp1 == 1) / (1)",  # extra parentheses
                ["mp1 == 1"],
                ["1"],
            ),
        ],
    )
    def test_parse(self, expr, numerators, denominators):
        name, num_list, den_list = fml._parse_formula_expression(expr)
        assert name == "A"
        assert [term.expr for term in num_list] == numerators
        assert [term.expr for term in den_list] == denominators

    @pytest.mark.parametrize(
        "expr",
        [
            "",  # Empty string
            "random string",  # Random string
            "A = (mp1 == 1) / (mp2 > 2) / 1",  # No more than one '/'
            "A = (mp1 == 1)",  # No denominator
            "(mp1 == 1) * 2 / 1",  # No name
            "A = mp1 == 1 * 2 / 1",  # No parentheses
            "A = (mp1 == 1) & (mp2 > 2) / 1",  # Invalid operator '&'
            "A = (mp1 = 1) / 1",  # Invalid operator '='
        ],
    )
    def test_invalid_expression(self, expr):
        with pytest.raises(fml.FormulaError):
            fml._parse_formula_expression(expr)


class TestTraitFormula:
    """Test `TraitFormula` class."""

    def test_init(self):
        description = "Some description"
        expression = "A = (mp1 == 1) / 1"
        formula = fml.TraitFormula(expression, description)
        assert formula.name == "A"
        assert formula.description == description
        assert formula.expression == expression
        assert [term.expr for term in formula._numerators] == ["mp1 == 1"]
        assert [term.expr for term in formula._denominators] == ["1"]

    def test_sia_linkage(self):
        description = "The ratio of sialylated to non-sialylated glycans"
        expression = "CS = nS // (type == 'complex')"
        formula = fml.TraitFormula(expression, description)
        assert formula.sia_linkage is True

        description = "The ratio of high-mannose to hybrid glycans"
        expression = "MHy = (type == 'high-mannose') / (type == 'hybrid')"
        formula = fml.TraitFormula(expression, description)
        assert formula.sia_linkage is False

    def test_numerators_and_denominators(self):
        description = "Some description"
        expression = "A = (mp1 == 1) / 1"
        formula = fml.TraitFormula(expression, description)
        assert formula.numerators == ["mp1 == 1"]
        assert formula.denominators == ["1"]

    @pytest.fixture
    def abund_table(self):
        #      G1  G2  G3
        # S1   1   2   2
        # S2   2   1   2
        # S3   1   2   1

        data = {
            "G1": [1, 2, 1],
            "G2": [2, 1, 2],
            "G3": [2, 2, 1],
        }
        return pd.DataFrame(data, index=["S1", "S2", "S3"])

    @pytest.mark.parametrize(
        "expression, expected",
        [
            (
                "A = (mp_int == 1) / 1",
                [0.2, 0.4, 0.25],
            ),
            (
                "A = (mp_int == 3) // (mp_int >= 2)",
                [0.5, 2 / 3, 1 / 3],
            ),
            (
                "A = mp_int // (mp_bool == True)",
                [7 / 3, 2, 2],
            ),
        ],
    )
    def test_calcu_trait(self, mp_table, abund_table, expression, expected):
        # mp_table:
        #     mp_int  mp_bool mp_str
        # G1       1     True      a
        # G2       2    False      b
        # G3       3     True      c
        formula = fml.TraitFormula(expression, "Some description")
        formula.initialize(mp_table)
        result = formula.calcu_trait(abund_table)
        expected = pd.Series(expected, index=abund_table.index, name="A")
        pd.testing.assert_series_equal(result, expected)


@pytest.mark.skip("`TraitFormula` to be updated.")
def test_load_default_formulas():
    structure_formulas = list(fml.load_default_formulas("structure"))
    composition_formulas = list(fml.load_default_formulas("composition"))
    assert len(structure_formulas) > 0
    assert len(composition_formulas) > 0
    assert structure_formulas[0].type == "structure"
    assert composition_formulas[0].type == "composition"
    assert len(structure_formulas) != len(composition_formulas)


class TestLoadFormulasFromFile:
    """Test `load_formulas_from_file` function."""

    def test_basic(self, write_content):
        description = "The ratio of high-mannose to hybrid glycans"
        expression = "MHy = (type == 'high-mannose') / (type == 'hybrid')"
        file = write_content(f"@ {description}\n$ {expression}")
        result = list(fml.load_formulas_from_file(file))
        assert len(result) == 1
        assert result[0].description == description
        assert result[0].name == "MHy"

    def test_duplicated_formulas(self, write_content):
        description1 = "The ratio of high-mannose to hybrid glycans"
        expression1 = "MHy = (type == 'high-mannose') / (type == 'hybrid')"
        description2 = (
            "Relative abundance of high mannose type glycans within total spectrum"
        )
        expression2 = (
            "MHy = (type == 'high-mannose') / (type == 'hybrid')"  # same expression
        )
        content = (
            f"@ {description1}\n$ {expression1}\n@ {description2}\n$ {expression2}"
        )
        file = write_content(content)
        with pytest.raises(FormulaError) as excinfo:
            list(fml.load_formulas_from_file(file))
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
