import numpy as np
import pandas as pd
import pytest
from attrs import define
from hypothesis import given, assume, strategies as st

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
            ("1", 1.0),
            ("(1)", 1.0),
            ("1.0", 1.0),
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

    @pytest.mark.parametrize("mp", ["mp_str"])
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
            ("mp_int > 2", "mp_int", ">", 2),
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


class TestDivisionTermWrapper:

    @given(
        values=st.lists(
            st.floats(
                allow_nan=False,
                allow_infinity=False,
                min_value=0.01,
                max_value=100.0,
            ),
            min_size=2,
        )
    )
    def test_call_values(self, values):
        @define
        class MockTerm:
            expr = "mock_term"

            def __call__(self, table):
                return pd.Series(
                    data=values,
                    name=self.expr,
                    dtype="Float32",
                    index=[f"G{i}" for i in range(len(values))],
                )

        term = MockTerm()
        wrapper = fml.DivisionTermWrapper(term)
        result = wrapper(None)
        assert np.allclose(result.values * values, 1.0, rtol=1e-6)

    @given(st.text(min_size=1))
    def test_call_expr(self, original_expr):
        @define
        class MockTerm:
            expr = original_expr
            __call__ = lambda self, table: pd.Series(
                [1, 2, 3], name=self.expr, dtype="Float32", index=["G1", "G2", "G3"]
            )

        term = MockTerm()
        wrapper = fml.DivisionTermWrapper(term)
        result = wrapper(None)
        assert result.name == f"/ ({term.expr})"

    def test_raise(self):
        @define
        class MockTerm:
            expr = "mock_term"

            def __call__(self, table):
                raise fml.FormulaTermCalculationError("error")

        term = MockTerm()
        wrapper = fml.DivisionTermWrapper(term)
        with pytest.raises(fml.FormulaTermCalculationError):
            wrapper(None)

    def test_zeros(self):
        @define
        class MockTerm:
            expr = "mock_term"

            def __call__(self, table):
                return pd.Series(
                    [1, 0, 3], index=["G1", "G2", "G3"], name=self.expr, dtype="Float32"
                )

        term = MockTerm()
        wrapper = fml.DivisionTermWrapper(term)
        result = wrapper(None)
        assert np.allclose(result.values, [1, 0, 1 / 3], rtol=1e-6)

    def test_from_compare_term(self):
        term = fml.CompareTerm("mp_int", ">", 2)
        with pytest.raises(ValueError):
            wrapper = fml.DivisionTermWrapper(term)


class TestParseFormulaExpression:
    """Test `_parse_formula_expression` function."""

    @given(name=st.text(min_size=1), num=st.text(min_size=1), den=st.text(min_size=1))
    @pytest.mark.parametrize("spliter", ["//", "/"])
    def test_split_expr(self, name, num, den, spliter):
        not_appeared = "/[]\n\r\t"
        assume(not any(c in name for c in not_appeared))
        assume(not any(c in num for c in not_appeared))
        assume(not any(c in den for c in not_appeared))
        assume(" " not in name)
        assume(name.strip() != "")
        assume(num.strip() != "")
        assume(den.strip() != "")

        expr = f"{name} = [{num}] {spliter} [{den}]"
        parser = fml.FormulaParser()
        result = parser._split_expr(expr)
        expected = (name.strip(), num.strip(), spliter, den.strip())
        assert result == expected

    @pytest.mark.parametrize(
        "expr, expected",
        [
            ("a * b", [("*", "a"), ("*", "b")]),
            ("a * b / c", [("*", "a"), ("*", "b"), ("/", "c")]),
            ("a / b / c", [("*", "a"), ("/", "b"), ("/", "c")]),
            ("a", [("*", "a")]),
            ("(a == 1) / b * c", [("*", "(a == 1)"), ("/", "b"), ("*", "c")]),
            ("a / 2", [("*", "a"), ("/", "2")]),
        ],
    )
    def test_split_terms(self, expr, expected):
        parser = fml.FormulaParser()
        result = parser._split_terms(expr)
        assert result == expected

    def test_parse_term(self):
        VALID_EXPR = "__valid_expr__"

        @define
        class FakeTerm:
            expr = "fake_term"

            @classmethod
            def from_expr(cls, expr):
                if expr == VALID_EXPR:
                    return cls()
                raise fml.FormulaTermParseError("Invalid expression")

        parser = fml.FormulaParser(available_terms=[FakeTerm])
        result = parser._parse_term(VALID_EXPR)
        assert result.expr == "fake_term"

    def test_parse_term_invalid(self):
        INVALID_EXPR = "__invalid_expr__"

        @define
        class FakeTerm:
            expr = "fake_term"

            @classmethod
            def from_expr(cls, expr):
                if expr == INVALID_EXPR:
                    raise fml.FormulaTermParseError("Invalid expression")
                return cls()

        parser = fml.FormulaParser(available_terms=[FakeTerm])
        with pytest.raises(fml.FormulaTermParseError):
            parser._parse_term(INVALID_EXPR)

    def test_parse_terms_expr(self, mocker):
        split_terms_result = [("*", "a"), ("/", "b")]
        mocker.patch(
            "glytrait.formula.FormulaParser._split_terms",
            return_value=split_terms_result,
        )

        term1 = mocker.Mock()
        term2 = mocker.Mock()
        terms = [term1, term2]
        mocker.patch(
            "glytrait.formula.FormulaParser._parse_term",
            side_effect=terms,
        )

        parser = fml.FormulaParser()
        r_term1, r_term2 = parser._parse_terms_expr("expr")
        assert r_term1 == term1
        assert r_term2 == fml.DivisionTermWrapper(term2)

    def test_parse_one_slash(self, mocker):
        mocker.patch(
            "glytrait.formula.FormulaParser._split_expr",
            return_value=("name", "num", "/", "den"),
        )
        mocker.patch(
            "glytrait.formula.FormulaParser._parse_terms_expr",
            side_effect=[["a", "b"], ["c", "d"]],
        )

        @define
        class FakeFormula:
            name: str
            numerators: list[str]
            denominators: list[str]

        parser = fml.FormulaParser(formula_factory=FakeFormula)
        result = parser.parse("expr")

        assert result.name == "name"
        assert result.numerators == ["a", "b"]
        assert result.denominators == ["c", "d"]

    def test_parse_two_slashes(self, mocker):
        mocker.patch(
            "glytrait.formula.FormulaParser._split_expr",
            return_value=("name", "num", "//", "den"),
        )
        mocker.patch(
            "glytrait.formula.FormulaParser._parse_terms_expr",
            side_effect=[["a", "b"], ["c", "d"]],
        )

        @define
        class FakeFormula:
            name: str
            numerators: list[str]
            denominators: list[str]

        parser = fml.FormulaParser(formula_factory=FakeFormula)
        result = parser.parse("expr")

        assert result.name == "name"
        assert result.numerators == ["a", "b", "c", "d"]
        assert result.denominators == ["c", "d"]


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


@pytest.mark.parametrize(
    "sia_linkage, expected", [(False, ["F1"]), (True, ["F1", "F2"])]
)
def test_load_formulas(mocker, sia_linkage, expected):
    @define
    class FakeFormula:
        name: str
        sia_linkage: bool

    @define
    class FakeParser:
        def parse(self, file):
            return [FakeFormula("F1", False), FakeFormula("F2", True)]

    mocker.patch("glytrait.formula.FormulaFileParser", return_value=FakeParser())
    formulas = fml.load_formulas("file", sia_linkage=sia_linkage)
    assert [f.name for f in formulas] == expected


def test_load_default_formulas():
    structure_formulas = list(fml.load_default_formulas("structure"))
    composition_formulas = list(fml.load_default_formulas("composition"))
    assert len(structure_formulas) > 0
    assert len(composition_formulas) > 0
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
