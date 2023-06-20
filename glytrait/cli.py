from pathlib import Path

import click
import emoji

from glytrait.config import Config
from glytrait.core import run_workflow
from glytrait.exception import GlyTraitError
from glytrait.trait import save_trait_formula_template

UNDIFINED = "__UNDEFINED__"


def save_template_callback(ctx, param, value):
    """Save a template for the user to fill in."""
    if value == UNDIFINED:
        return
    if Path(value).exists() and not Path(value).is_dir():
        raise click.BadParameter("The path to save the template must be a directory.")
    else:
        save_trait_formula_template(value)
        msg = (
            f"Template saved to {value}. :victory_hand:\n"
            f"You can edit the template and use `glytrait INPUT_FILE OUTPUT_FILE -f "
            f"TEMPLATE_FILE` to provide additional traits to glyTrait."
        )
        click.echo(emoji.emojize(msg))
        ctx.exit()


def file_validator(suffix: str, file_name: str):
    """Validate file suffix."""

    def validator(ctx, param, value):
        if value is not None and Path(value).suffix != suffix:
            raise click.BadParameter(
                f"{file_name.capitalize()} must be a {suffix} file."
            )
        return value

    return validator


@click.command()
@click.argument(
    "input_file",
    type=click.Path(exists=True),
    required=False,
    callback=file_validator(".csv", "input file"),
)
@click.option(
    "-o",
    "--output_file",
    type=click.Path(),
    callback=file_validator(".xlsx", "output file"),
)
@click.option(
    "-r",
    "--filter_glycan_ratio",
    type=click.FLOAT,
    default=0.5,
    help="Glycans with missing value ratio greater than this value will be filtered out.",
)
@click.option(
    "-i",
    "--impute_method",
    type=click.Choice(["zero", "min", "lod", "mean", "median"]),
    default="min",
    help="Method to impute missing values.",
)
@click.option(
    "--sia_linkage", is_flag=True, help="Include sialic acid linkage traits."
)
@click.option(
    "-f",
    "--formula_file",
    type=click.Path(exists=True),
    help="User formula file.",
    callback=file_validator(".txt", "formula file"),
)
@click.option(
    "-t",
    "--save_template",
    type=click.Path(),
    callback=save_template_callback,
    is_eager=True,
    expose_value=False,
    default=UNDIFINED,
    help="Save a template for the user to fill in.",
)
@click.option(
    "--filter/--no-filter",
    default=True,
    help="Filter out invalid derived traits.",
)
@click.option(
    "-g",
    "--group_file",
    type=click.Path(exists=True),
    help="Group file for hypothesis test.",
    callback=file_validator(".csv", "group file"),
)
@click.option(
    "-s",
    "--structure_file",
    type=click.Path(exists=True),
    help="Structure file for hypothesis test.",
)
@click.option(
    "-d",
    "--database",
    type=click.STRING,
    help="Built-in database to use, either 'serum' or 'IgG'."
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
def cli(
    input_file,
    output_file,
    filter_glycan_ratio,
    impute_method,
    sia_linkage,
    formula_file,
    filter,
    group_file,
    structure_file,
    database,
    mode
):
    """Run the glytrait workflow."""
    if output_file is None:
        output_file = str(
            Path(input_file).with_name(Path(input_file).stem + "_glytrait.xlsx")
        )
    else:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    mode = "composition" if mode.lower() in ["c", "composition"] else "structure"
    try:
        config = Config(dict(
            input_file=input_file,
            output_file=output_file,
            mode=mode,
            filter_glycan_max_na=filter_glycan_ratio,
            impute_method=impute_method,
            sia_linkage=sia_linkage,
            formula_file=formula_file,
            filter_invalid_traits=filter,
            group_file=group_file,
            structure_file=structure_file,
            database=database,
        ))
        run_workflow(config)
    except GlyTraitError as e:
        raise click.UsageError(str(e) + emoji.emojize(" :thumbs_down:"))
    msg = f"Done :thumbs_up:! Output written to {output_file}."
    click.echo(emoji.emojize(msg))


if __name__ == "__main__":
    cli()
