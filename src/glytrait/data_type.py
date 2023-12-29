from __future__ import annotations

from typing import NewType

import pandas as pd

AbundanceTable = NewType("AbundanceTable", pd.DataFrame)
"""Abundance table type.
The index are samples and the columns are glycans.
The abundance table could only be returned by `load_abundance_table` function.
"""

GroupSeries = NewType("GroupSeries", pd.Series)
"""Group series type.
The index are the sample names and the values are groups.
"""

MetaPropertyTable = NewType("MetaPropertyTable", pd.DataFrame)
"""The type of the meta-property table. 
Only returned by `build_meta_property_table`.
A pandas DataFrame with `glycan_ids` as the index, 
and the meta-property names as the columns.
"""

DerivedTraitTable = NewType("DerivedTraitTable", pd.DataFrame)
"""The type of the derived trait table.
Only returned by `calcu_derived_trait`.
A pandas DataFrame with samples as the index,
and the derived trait names as the columns.
"""
