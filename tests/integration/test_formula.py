import pytest

import glytrait.formula as fml


@pytest.mark.parametrize(
    "expr, numerators, denominators",
    [
        (
            "T = [A] / [B]",
            {"A"},
            {"B"},
        ),
        (
            "T = [A] // [B]",
            {"A", "B"},
            {"B"},
        ),
        (
            "T = [A * B] / [C]",
            {"A", "B"},
            {"C"},
        ),
        (
            "T = [A * B] / [C * D]",
            {"A", "B"},
            {"C", "D"},
        ),
        (
            "T = [A / B] / [C]",
            {"A", "/ (B)"},
            {"C"},
        ),
        (
            "T = [(A == 1)] // [B]",
            {"A == 1", "B"},
            {"B"},
        ),
        (
            "T = [(A != 1)] // [B]",
            {"A != 1", "B"},
            {"B"},
        ),
        (
            "T = [(A > 1)] // [B]",
            {"A > 1", "B"},
            {"B"},
        ),
        (
            "T = [(A >= 1)] // [B]",
            {"A >= 1", "B"},
            {"B"},
        ),
        (
            "T = [(A < 1)] // [B]",
            {"A < 1", "B"},
            {"B"},
        ),
        (
            "T = [(A <= 1)] // [B]",
            {"A <= 1", "B"},
            {"B"},
        ),
        (
            "T = [(A == 1) * (B == 2)] // [C]",
            {"A == 1", "B == 2", "C"},
            {"C"},
        ),
        (
            "T = [A / B / C] / [D]",
            {"A", "/ (B)", "/ (C)"},
            {"D"},
        ),
        (
            "T = [A / B * C] / [D]",
            {"A", "/ (B)", "C"},
            {"D"},
        ),
        (
            "T = [A] / [1]",
            {"A"},
            {"1.0"},
        ),
        (
            "T = [A] / [1.0]",
            {"A"},
            {"1.0"},
        ),
        (
            "T = [A / 2] / [1]",
            {"A", "/ (2.0)"},
            {"1.0"},
        ),
    ],
)
def test_parse_formula(expr, numerators, denominators):
    parser = fml.FormulaParser()
    result_f = parser.parse(expr)
    assert result_f.name == "T"
    assert set(t.expr for t in result_f.numerators) == numerators
    assert set(t.expr for t in result_f.denominators) == denominators


def test_load_default_formulas():
    structure_formulas = list(fml.load_default_formulas("structure"))
    composition_formulas = list(fml.load_default_formulas("composition"))
    assert len(structure_formulas) > 0
    assert len(composition_formulas) > 0
    assert len(structure_formulas) != len(composition_formulas)
