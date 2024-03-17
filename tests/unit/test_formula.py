import pandas as pd
import pytest

import glytrait.formula as fml


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


def mock_formula_parser(expr: str):
    """Mock formula parser."""
    return


class TestConstantTerm:

    def test_expr(self):
        term = fml.ConstantTerm(1)
        assert term.expr == "1"

    def test_call(self, mp_table):
        term = fml.ConstantTerm(1)
        result = term(mp_table)
        expected = pd.Series([1, 1, 1], index=mp_table.index, name="1", dtype="Float32")
        pd.testing.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "expr, value",
        [
            ("1", 1.),
            ("(1)", 1.),
            ("1.0", 1.),
            ("0.5", 0.5),
        ],
    )
    def test_from_expr(self, expr, value):
        term = fml.ConstantTerm.from_expr(expr)
        assert term.value == value

    @pytest.mark.parametrize("expr", ["a"])
    def test_from_expr_invalid(self, expr):
        with pytest.raises(fml.FormulaParseError):
            fml.ConstantTerm.from_expr(expr)

    @pytest.mark.parametrize("value", [0, -1])
    def test_value_error(self, value):
        with pytest.raises(ValueError):
            fml.ConstantTerm(value)


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
        with pytest.raises(fml.MetaPropertyTypeError):
            term(mp_table)

    def test_call_mp_not_exist(self, mp_table):
        term = fml.NumericalTerm("mp_not_exist")
        with pytest.raises(fml.MissingMetaPropertyError):
            term(mp_table)

    @pytest.mark.parametrize(
        "expr, mp",
        [
            ("mp", "mp"),
            ("(mp)", "mp"),
        ],
    )
    def test_from_expr(self, expr, mp):
        term = fml.NumericalTerm.from_expr(expr)
        assert term.meta_property == mp

    @pytest.mark.parametrize("expr", ["1", "mp > 2", "mp == 2"])
    def test_from_expr_invalid(self, expr):
        with pytest.raises(fml.FormulaParseError):
            fml.NumericalTerm.from_expr(expr)


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
        with pytest.raises(fml.MetaPropertyTypeError):
            term(mp_table)

    @pytest.mark.parametrize(
        "expr, mp, operator, value",
        [
            ("(mp_int > 2)", "mp_int", ">", 2),
            ("(mp_int >= 2)", "mp_int", ">=", 2),
            ("(mp_int < 2)", "mp_int", "<", 2),
            ("(mp_int <= 2)", "mp_int", "<=", 2),
            ("(mp_int == 2)", "mp_int", "==", 2),
            ("(mp_int != 2)", "mp_int", "!=", 2),
            ("(mp_int > 2)", "mp_int", ">", 2),
            ("(mp_bool == True)", "mp_bool", "==", True),
            ("(mp_str != 'b')", "mp_str", "!=", "b"),
            ('(mp_str != "b")', "mp_str", "!=", "b"),
        ],
    )
    def test_from_expr(self, expr, mp, operator, value):
        term = fml.CompareTerm.from_expr(expr)
        assert term.meta_property == mp
        assert term.operator == operator
        assert term.value == value

    @pytest.mark.parametrize(
        "expr",
        [
            "mp_int > 2",  # no parentheses
            "(mp_int > 2) > 1",  # more than one '>'
            "(mp_int > )",  # missing value
            "(mp_bool == True) & (mp_str != 'b')",  # invalid operator '&'
            "(mp_str)",  # no operator
            "(1)",  # invalid meta property
        ],
    )
    def test_from_expr_invalid(self, expr):
        with pytest.raises(fml.FormulaTermParseError):
            fml.CompareTerm.from_expr(expr)

    def test_missing_mp(self, mp_table):
        term = fml.CompareTerm("mp_not_exist", ">", 2)
        with pytest.raises(fml.MissingMetaPropertyError):
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
                ["1.0"],
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
                ["1.0"],
            ),
            (
                "A = (mp1 == 1) * (mp2 > 2) * (mp3 <= 4) / 1",
                ["mp1 == 1", "mp2 > 2", "mp3 <= 4"],
                ["1.0"],
            ),
            (
                "A = (mp1 == 1) * 2 / 1",
                ["mp1 == 1", "2.0"],
                ["1.0"],
            ),
            (
                "A = mp1 / mp2",
                ["mp1"],
                ["mp2"],
            ),
            (
                "A = (mp1 == 1) / (1)",  # extra parentheses
                ["mp1 == 1"],
                ["1.0"],
            ),
        ],
    )
    def test_parse(self, expr, numerators, denominators):
        parser = fml.FormulaParser()
        formula = parser(expr)
        assert formula.name == "A"
        assert [term.expr for term in formula.numerators] == numerators
        assert [term.expr for term in formula.denominators] == denominators

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
        parser = fml.FormulaParser()
        with pytest.raises(fml.FormulaParseError):
            parser(expr)


class TestTraitFormula:
    """Test `TraitFormula` class."""

    @pytest.mark.parametrize(
        "numerators, denominators, expected",
        [
            (
                [fml.NumericalTerm("nS")],
                [fml.NumericalTerm("1")],
                False,
            ),
            (
                [fml.NumericalTerm("nL")],
                [fml.NumericalTerm("1")],
                True,
            ),
            (
                [fml.NumericalTerm("nE")],
                [fml.NumericalTerm("1")],
                True,
            ),
        ],
        ids=["nS", "nL", "nE"],
    )
    def test_sia_linkage(self, numerators, denominators, expected):
        formula = fml.TraitFormula("F", numerators, denominators)
        assert formula.sia_linkage == expected

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

    def test_initialize_failed(self, mp_table):
        numerators = [fml.NumericalTerm("mp_not_exist")]
        denominators = [fml.NumericalTerm("1")]
        formula = fml.TraitFormula("F", numerators, denominators)
        with pytest.raises(fml.FormulaCalculationError):
            formula.initialize(mp_table)

    def test_calcu_trait_without_initialization(self, mp_table, abund_table):
        formula = fml.TraitFormula(
            name="F",
            numerators=[fml.NumericalTerm("mp_int")],
            denominators=[fml.NumericalTerm("1")],
        )
        with pytest.raises(fml.FormulaNotInitializedError):
            formula.calcu_trait(abund_table)

    @pytest.mark.parametrize(
        "numerators, denominators, expected",
        [
            (
                [fml.CompareTerm("mp_int", "==", 1)],
                [fml.ConstantTerm(1)],
                [0.2, 0.4, 0.25],
            ),
            (
                [
                    fml.CompareTerm("mp_int", "==", 3),
                    fml.CompareTerm("mp_int", ">=", 2),
                ],
                [fml.CompareTerm("mp_int", ">=", 2)],
                [0.5, 2 / 3, 1 / 3],
            ),
            (
                [fml.NumericalTerm("mp_int"), fml.CompareTerm("mp_bool", "==", True)],
                [fml.CompareTerm("mp_bool", "==", True)],
                [7 / 3, 2, 2],
            ),
        ],
    )
    def test_calcu_trait(
        self, mp_table, abund_table, numerators, denominators, expected
    ):
        # mp_table:
        #     mp_int  mp_bool mp_str
        # G1       1     True      a
        # G2       2    False      b
        # G3       3     True      c
        formula = fml.TraitFormula("F", numerators, denominators)
        formula.initialize(mp_table)
        result = formula.calcu_trait(abund_table)
        expected = pd.Series(expected, index=abund_table.index, name="F")
        pd.testing.assert_series_equal(result, expected)

    def test_from_expr(self, mocker):
        parser = mocker.Mock()
        parser.return_value = "formula"
        result = fml.TraitFormula.from_expr("expr", parser=parser)
        assert result == "formula"
        parser.assert_called_once_with("expr")


@pytest.mark.skip("`TraitFormula` to be updated.")
def test_load_default_formulas():
    structure_formulas = list(fml.load_default_formulas("structure"))
    composition_formulas = list(fml.load_default_formulas("composition"))
    assert len(structure_formulas) > 0
    assert len(composition_formulas) > 0
    assert structure_formulas[0].type == "structure"
    assert composition_formulas[0].type == "composition"
    assert len(structure_formulas) != len(composition_formulas)


def test_save_builtin_formula(clean_dir):
    fml.save_builtin_formula(clean_dir)
    struc_file = clean_dir / "struc_builtin_formulas.txt"
    comp_file = clean_dir / "comp_builtin_formulas.txt"
    assert struc_file.exists()
    assert comp_file.exists()


class TestFormulaFileParser:

    def test_parse(self, write_content):
        file = write_content("@ Description\n$ Expression\n")
        parser = fml.FormulaFileParser(expr_parser=lambda x: x)
        assert list(parser.parse(file)) == ["Expression"]

    def test_deconvolute_formula_file_basic(self, write_content):
        file = write_content("@ Description\n$ Expression\n")
        result = list(fml.FormulaFileParser._deconvolute_formula_file(file))
        expected = [("Description", "Expression")]
        assert result == expected

    def test_first_line_not_description(self, write_content):
        file = write_content("$ Expression\n@ Description\n")
        with pytest.raises(fml.FormulaFileError) as excinfo:
            list(fml.FormulaFileParser._deconvolute_formula_file(file))
        assert "No description before expression 'Expression'" in str(excinfo.value)

    def test_two_descriptions(self, write_content):
        file = write_content("@ Description1\n@ Description2\n")
        with pytest.raises(fml.FormulaFileError) as excinfo:
            list(fml.FormulaFileParser._deconvolute_formula_file(file))
        assert "No expression follows description 'Description1'." in str(excinfo.value)

    def test_two_expressions(self, write_content):
        file = write_content("@ Description\n$ Expression1\n$ Expression2")
        with pytest.raises(fml.FormulaFileError) as excinfo:
            list(fml.FormulaFileParser._deconvolute_formula_file(file))
        assert "No description before expression 'Expression2'." in str(excinfo.value)

    def test_no_last_expression(self, write_content):
        file = write_content("@ Description1\n$ Expression1\n@Description2")
        with pytest.raises(fml.FormulaFileError) as excinfo:
            list(fml.FormulaFileParser._deconvolute_formula_file(file))
        assert "No expression follows description 'Description2'." in str(excinfo.value)
