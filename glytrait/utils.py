from glypy.structure.glycan_composition import MonosaccharideResidue


def get_mono_comp(mono: MonosaccharideResidue) -> str:
    """Get the composition of a monosaccharide residue.

    Args:
        mono (MonosaccharideResidue): The monosaccharide residue.

    Returns:
        GlycanComposition: The composition.
    """
    return MonosaccharideResidue.from_monosaccharide(mono).name()
