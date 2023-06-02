from pathlib import Path

import click
import emoji

from glytrait.core import run_workflow
from glytrait.exception import GlyTraitError
from glytrait.trait import (
    save_trait_formula_template
)

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


def valid_input_file(ctx, param, value):
    """Validate input file."""
    if Path(value).suffix != ".csv":
        raise click.BadParameter("Input file must be a csv file.")
    return value


def valid_output_file(ctx, param, value):
    """Validate output file."""
    if value is not None and Path(value).suffix != ".xlsx":
        raise click.BadParameter("Output file must be a xlsx file.")
    return value


def valid_formula_file(ctx, param, value):
    """Validate formula file."""
    if value is not None and Path(value).suffix != ".txt":
        raise click.BadParameter("Formula file must be a txt file.")
    return value


@click.command()
@click.argument(
    "input_file",
    type=click.Path(exists=True),
    required=False,
    callback=valid_input_file,
)
@click.option("-o", "--output_file", type=click.Path(), callback=valid_output_file)
@click.option(
    "-s", "--sia_linkage", is_flag=True, help="Include sialic acid linkage traits."
)
@click.option(
    "-f",
    "--formula_file",
    type=click.Path(exists=True),
    help="User formula file.",
    callback=valid_formula_file,
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
def cli(input_file, output_file, sia_linkage, formula_file):
    """Run the glytrait workflow.

    You can use the `--save_template` option to save the formula template
    to the specified directory, then edit the template and
    use it to provide additional traits to glyTrait.
    """
    if output_file is None:
        output_file = str(
            Path(input_file).with_name(Path(input_file).stem + "_glytrait.xlsx")
        )
    else:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    try:
        run_workflow(input_file, output_file, sia_linkage, formula_file)
    except GlyTraitError as e:
        raise click.UsageError(str(e) + emoji.emojize(" :thumbs_down:"))
    msg = f"Done :thumbs_up:! Output written to {output_file}."
    click.echo(emoji.emojize(msg))


if __name__ == "__main__":
    cli()
