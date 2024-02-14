import pandas as pd
import pytest

from glytrait.data_export import export_all
from glytrait.formula import TraitFormula


class TestExportAll:
    def test_export_dataframe(self, clean_dir):
        df = pd.DataFrame(
            {
                "a": [1, 2],
                "b": [3, 4],
            },
            index=pd.Index([1, 2], name="index"),
        )
        filename = "test.csv"

        export_all([(filename, df)], clean_dir)

        assert (clean_dir / filename).exists()
        assert (clean_dir / filename).read_text() == ("index,a,b\n" "1,1,3\n" "2,2,4\n")

    @pytest.mark.skip("`TraitFormula` to be updated.")
    def test_export_formulas(self, clean_dir):
        formula_1 = TraitFormula(
            description="The ratio of high-mannose to complex glycans",
            name="MHy",
            type="structure",
            numerator_properties=["isHighMannose"],
            denominator_properties=["isComplex"],
        )
        formula_2 = TraitFormula(
            description="The ratio of complex glycans",
            name="TC",
            type="structure",
            numerator_properties=["isComplex"],
            denominator_properties=["."],
        )
        formulas = [formula_1, formula_2]
        filename = "test.csv"

        export_all([(filename, formulas)], clean_dir)

        assert (clean_dir / filename).exists()
        assert (clean_dir / filename).read_text() == (
            "MHy: The ratio of high-mannose to complex glycans\n"
            "TC: The ratio of complex glycans\n"
        )
