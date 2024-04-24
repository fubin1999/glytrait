"""Command line interface for glyTrait."""

from pathlib import Path
from typing import Literal, Any

import click
import emoji
import pandas as pd

from glytrait.load_data import load_data
from glytrait.api import Experiment, MissingDataError
from glytrait.formula import save_builtin_formula
from glytrait.data_export import export_all

UNDIFINED = "__UNDEFINED__"


def save_builtin_formulas_callback(ctx, param, value):
    """Save a copy of the built-in formulas."""
    if value == UNDIFINED:
        return
    if Path(value).exists() and not Path(value).is_dir():
        msg = "The path to save the built-in formulas must be a directory."
        raise click.BadParameter(msg)
    else:
        save_builtin_formula(value)
        msg = (
            f"Built-in formulas saved to: "
            f"{value}/struc_builtin_formulas.txt, {value}/comp_builtin_formulas.txt"
        )
        click.echo(emoji.emojize(msg))
        ctx.exit()


@click.command()
@click.argument(
    "abundance-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=False,
)
@click.argument(
    "glycan-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
    required=False,
)
@click.option(
    "-g",
    "--group-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="Group file path.",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(dir_okay=True, file_okay=False),
    help="Output file path. Default is the input file name with '_glytrait.xlsx' "
    "suffix.",
)
@click.option(
    "-m",
    "--mode",
    type=click.Choice(["structure", "composition", "S", "C"]),
    default="structure",
    help="Mode to run glyTrait, either 'structure' or 'composition'. "
    "You can also use 'S' or 'C' for short. "
    "Default is 'structure'.",
)
@click.option(
    "-r",
    "--filter-glycan-ratio",
    type=click.FLOAT,
    default=1,
    help="Glycans with missing value ratio greater than this value will be filtered out.",
)
@click.option(
    "-i",
    "--impute-method",
    type=click.Choice(["zero", "min", "lod", "mean", "median"]),
    default="zero",
    help="Method to impute missing values.",
)
@click.option(
    "-l",
    "--sia-linkage",
    is_flag=True,
    help="Include sialic acid linkage traits.",
)
@click.option(
    "-f",
    "--formula-file",
    type=click.Path(exists=True),
    help="User formula file.",
)
@click.option(
    "--filter/--no-filter",
    default=True,
    help="Filter out invalid derived traits. Default is filtering."
    "Use --no-filter to disable filtering.",
)
@click.option(
    "-c",
    "--corr-threshold",
    type=click.FLOAT,
    default=1,
    help="Threshold for correlation between traits. "
    "Default is 1, which means only traits with perfect collinearity "
    "will be filtered. Use -1 to disable filtering.",
)
@click.option(
    "-b",
    "--builtin-formulas",
    type=click.Path(),
    callback=save_builtin_formulas_callback,
    is_eager=True,
    expose_value=False,
    default=UNDIFINED,
    help="The directory path to save a copy of the built-in formulas.",
)
@click.version_option()
def cli(
    abundance_file,
    glycan_file,
    group_file,
    output_dir,
    mode,
    filter_glycan_ratio,
    impute_method,
    sia_linkage,
    formula_file,
    filter,
    corr_threshold,
):
    """Run the glytrait workflow."""
    if abundance_file is None:
        _show_welcome_msg()
        return

    if output_dir is None:
        output_dir = _default_output_dir(abundance_file)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    mode = "composition" if mode.lower() in ["c", "composition"] else "structure"
    input_data = _make_input_data(abundance_file, glycan_file, group_file, mode)
    exp = Experiment(input_data, sia_linkage=sia_linkage, mode=mode)

    exp.preprocess(filter_glycan_ratio, impute_method)
    exp.extract_meta_properties()
    exp.derive_traits()
    if filter:
        exp.post_filter(corr_threshold)
        if group_file:
            exp.diff_analysis()

    output_data = _prepare_output(exp)
    export_all(output_data, output_dir)

    success_msg = f"Done :thumbs_up:! Output written to {output_dir}."
    click.echo(emoji.emojize(success_msg))


def _show_welcome_msg() -> None:
    msg = r"""
    Welcome to GlyTrait!

       _____ _    _______        _ _   
      / ____| |  |__   __|      (_) |  
     | |  __| |_   _| |_ __ __ _ _| |_ 
     | | |_ | | | | | | '__/ _` | | __|
     | |__| | | |_| | | | | (_| | | |_ 
      \_____|_|\__, |_|_|  \__,_|_|\__|
                __/ |                  
               |___/                   

    Use `glytrait --help` for more information.
    """
    click.echo(msg)


def _default_output_dir(abundance_file) -> str:
    dirpath = Path(abundance_file).with_name(Path(abundance_file).stem + "_glytrait")
    return str(dirpath)


def _make_input_data(
    abundance_file: str,
    glycan_file: str,
    group_file: str | None,
    mode: Literal["structure", "composition"],
):
    abundance_df = pd.read_csv(abundance_file)
    glycan_df = pd.read_csv(glycan_file)
    group_df = pd.read_csv(group_file) if group_file else None
    return load_data(abundance_df, glycan_df, group_df, mode=mode)


def _prepare_output(exp: Experiment) -> list[tuple[str, Any]]:
    result = {
        "derived_traits.csv": exp.derived_trait_table,
        "derived_traits_filtered.csv": exp.filtered_derived_trait_table,
        "glycan_abundance_processed.csv": exp.processed_abundance_table,
        "meta_properties.csv": exp.meta_property_table,
    }
    try:
        diff_results = exp.diff_results
    except MissingDataError:
        pass
    else:
        if "t_test" in diff_results:
            result["t_test.csv"] = diff_results["t_test"]
        else:  # anova
            result["anova.csv"] = diff_results["anova"]
            result["post_hoc.csv"] = diff_results["post_hoc"]
    return list(result.items())


if __name__ == "__main__":
    cli()
