import pandas as pd

from glytrait import Experiment

from .. import glycoct as ct

# ========== H5N4S2 ==========
# Neu5Ac - Gal - Glc2NAc - Man
#                             \
#                              Man - Glc2NAc - Glc2NAc
#                             /
# Neu5Ac - Gal - Glc2NAc - Man

# ========== H5N5F1S1 ==========
# Neu5Ac - Gal - Glc2NAc - Man                   Fuc
#                             \                   |
#                    Glc2NAc - Man - Glc2NAc - Glc2NAc
#                             /
#          Gal - Glc2NAc - Man

# ========== H5N2 ==========
# Man
#    \
#     Man
#    /   \
# Man     Man - Glc2NAc - Glc2NAc
#        /
#     Man

# ========== H6N6F1S1 ==========
# Neu5Ac - Gal - Glc2NAc
#                       \
#                        Man
#                       /   \
#          Gal - Glc2NAc     \                   Fuc
#                             \                   |
#                    Glc2NAc - Man - Glc2NAc - Glc2NAc
#                             /
#          Gal - Glc2NAc - Man

# ========== H4N4F3S1 ==========
#                  Fuc
#                   |
# Neu5Ac - Gal - Glc2NAc - Man                   Fuc
#                             \                   |
#                  Fuc         Man - Glc2NAc - Glc2NAc
#                   |         /
#                Glc2NAc - Man

GLYCANS = [
    ("H5N4S2", ct.test_glycoct_1),
    ("H5N5F1S1", ct.test_glycoct_2),
    ("H5N2", ct.test_glycoct_3),
    ("H6N6F1S1", ct.test_glycoct_8),
    ("H4N4F3S1", ct.test_glycoct_9),
]
glycan_df = pd.DataFrame.from_records(GLYCANS, columns=["GlycanID", "Structure"])

abundance_df = pd.DataFrame(
    {
        "Sample": ["S1", "S2", "S3", "S4", "S5"],
        "H5N4S2": [1.0, 2.0, 3.0, 4.0, 5.0],
        "H5N5F1S1": [2.0, 4.0, 1.0, 6.0, 2.0],
        "H5N2": [5.0, 4.0, 3.0, 2.0, 1.0],
        "H6N6F1S1": [4.0, 1.0, 2.0, 1.0, 5.0],
        "H4N4F3S1": [1.0, 2.0, 1.0, 2.0, 1.0],
    }
)


def test_calculate_traits(tmp_path):
    abundance_df.to_csv(tmp_path / "abundance.csv", index=False)
    glycan_df.to_csv(tmp_path / "glycans.csv", index=False)
    experiment = Experiment(
        abundance_file=str(tmp_path / "abundance.csv"),
        glycan_file=str(tmp_path / "glycans.csv"),
    )
    experiment.preprocess()
    experiment.extract_meta_properties()

    expressions = [
        "CS = [nS] // [type == 'complex']",
        "CGS = [nS / nG] // [type == 'complex']",
        "A2Fc = [nFc > 0] // [(type == 'complex') * (nAnt == 2)]",
    ]
    result = experiment.try_formulas(expressions, squeeze=False)

    ad = abundance_df
    expected = pd.DataFrame(
        {
            "CS": (
                (ad["H5N4S2"] * 2 + ad["H5N5F1S1"] + ad["H6N6F1S1"] + ad["H4N4F3S1"])
                / (ad["H5N4S2"] + ad["H5N5F1S1"] + ad["H6N6F1S1"] + ad["H4N4F3S1"])
            ).values,
            "CGS": (
                (
                    ad["H5N4S2"]
                    + ad["H5N5F1S1"] / 2
                    + ad["H6N6F1S1"] / 3
                    + ad["H4N4F3S1"]
                )
                / (ad["H5N4S2"] + ad["H5N5F1S1"] + ad["H6N6F1S1"] + ad["H4N4F3S1"])
            ).values,
            "A2Fc": (
                (ad["H5N5F1S1"] + ad["H4N4F3S1"])
                / (ad["H5N5F1S1"] + ad["H4N4F3S1"] + ad["H5N4S2"])
            ).values,
        },
        index=pd.Index(["S1", "S2", "S3", "S4", "S5"], name="Sample"),
    )
    pd.testing.assert_frame_equal(result, expected)
