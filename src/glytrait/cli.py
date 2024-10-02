"""Command line interface for glyTrait."""

from pathlib import Path
from typing import Literal, Any

import click
import emoji
from attrs import define

from glytrait.exception import GlyTraitError
from glytrait.api import Experiment, MissingDataError
from glytrait.formula import save_builtin_formula, load_formulas_from_file
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
@click.option(
    "--glycan-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="Glycan structure or composition file.",
)
@click.option(
    "--mp-file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="Meta-property file.",
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
    help="Output directory path. Default is the input file name with '_glytrait' "
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
    mp_file,
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

    _warn_about_filtering(filter, group_file)
    output_dir = _set_output_dir(output_dir, abundance_file)
    mode = _determine_mode(mode)

    config = WorkflowConfig(
        mode=mode,
        sia_linkage=sia_linkage,
        filter_glycan_ratio=filter_glycan_ratio,
        impute_method=impute_method,
        filter=filter,
        corr_threshold=corr_threshold,
    )
    try:
        _run_workflow(
            abundance_file,
            glycan_file,
            mp_file,
            group_file,
            formula_file,
            config,
            output_dir,
        )
    except GlyTraitError as e:
        click.echo(f"Error: {e}")
        return
    success_msg = f"Done :thumbs_up:! Output written to {output_dir}."
    click.echo(emoji.emojize(success_msg))


@define(kw_only=True)
class WorkflowConfig:
    """Configuration for the glyTrait workflow."""

    mode: Literal["structure", "composition"]
    sia_linkage: bool
    filter_glycan_ratio: float
    impute_method: Literal["zero", "min", "lod", "mean", "median"]
    filter: bool
    corr_threshold: float


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


def _warn_about_filtering(filter: bool, group_file: str) -> None:
    if not filter and group_file is not None:
        click.echo(
            "Warning: differential analysis will be disabled when --no-filter is used."
        )


def _set_output_dir(output_dir: str | None, abundance_file: str) -> str:
    if output_dir is None:
        output_dir_path = Path(abundance_file).with_name(
            Path(abundance_file).stem + "_glytrait"
        )
    else:
        output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    return str(output_dir_path)


def _determine_mode(mode: str) -> Literal["structure", "composition"]:
    return "composition" if mode.lower() in ["c", "composition"] else "structure"


def _prepare_output(exp: Experiment) -> list[tuple[str, Any]]:
    result = {
        "derived_traits.csv": exp.derived_trait_table,
        "glycan_abundance_processed.csv": exp.processed_abundance_table,
        "meta_properties.csv": exp.meta_property_table,
    }

    try:
        result["derived_traits_filtered.csv"] = exp.filtered_derived_trait_table
    except MissingDataError:
        pass

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


def _run_workflow(
    abundance_file: str,
    glycan_file: str | None,
    mp_file: str | None,
    group_file: str | None,
    formula_file: str | None,
    config: WorkflowConfig,
    output_dir: str,
):
    exp = Experiment(
        abundance_file=abundance_file,
        glycan_file=glycan_file,
        meta_property_file=mp_file,
        group_file=group_file,
        sia_linkage=config.sia_linkage,
        mode=config.mode,
    )
    _process_exp(exp, formula_file, config)
    output_data = _prepare_output(exp)
    export_all(output_data, output_dir)


def _process_exp(
    exp: Experiment,
    formula_file: str | None,
    config: WorkflowConfig,
) -> None:
    exp.preprocess(config.filter_glycan_ratio, config.impute_method)

    if formula_file:
        formulas = load_formulas_from_file(formula_file, config.sia_linkage)
        exp.derive_traits(formulas)
    else:
        exp.derive_traits()

    if config.filter:
        exp.post_filter(config.corr_threshold)
        if exp.groups is not None:
            exp.diff_analysis()


if __name__ == "__main__":
    cli()
