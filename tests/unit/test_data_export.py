import pandas as pd

from glytrait.data_export import export_all


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
