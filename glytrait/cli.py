from typing import Optional

import click
import emoji

from glytrait.io import read_input, write_output
from glytrait.trait import (
    save_trait_formula_template,
    load_formulas,
    build_meta_property_table,
    calcu_trait,
)

UNDIFINED = "__UNDEFINED__"


def save_template_callback(ctx, param, value):
    """Save a template for the user to fill in."""
    if value == UNDIFINED:
        return
    else:
        save_trait_formula_template(value)
        msg = (
            f"Template saved to {value}. :victory_hand:\n"
            f"You can edit the template and use `glytrait input_file output_file -f "
            f"template_file` to provide additional traits to glyTrait."
        )
        click.echo(emoji.emojize(msg))
        ctx.exit()


@click.command()
@click.argument("input_file", type=click.Path(exists=True), required=False)
@click.argument("output_file", type=click.Path(), required=False)
@click.option(
    "-s", "--sia_linkage", is_flag=True, help="Include sialic acid linkage traits."
)
@click.option(
    "-f", "--formula_file", type=click.Path(exists=True), help="User formula file."
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

    You can use `--save_template` option to save the formula template
    to the specified directory. You can then edit the template and
    use it to provide additional traits to glyTrait.
    """
    if input_file is None or output_file is None:
        raise click.UsageError(
            "You must provide both an input file and an output file."
        )
    run_workflow(input_file, output_file, sia_linkage, formula_file)
    msg = f"Done :thumbs_up:! Output written to {output_file}."
    click.echo(emoji.emojize(msg))


def run_workflow(
    input_file: str,
    output_file: str,
    sia_linkage: bool = False,
    user_formula_file: Optional[str] = None,
) -> None:
    """Run the workflow."""
    glycans, abund_df = read_input(input_file)
    formulas = load_formulas(user_formula_file)
    if not sia_linkage:
        formulas = [f for f in formulas if f.sia_linkage is False]
    meta_prop_df = build_meta_property_table(abund_df.columns, glycans, sia_linkage)
    trait_df = calcu_trait(abund_df, meta_prop_df, formulas)
    write_output(output_file, trait_df, abund_df, meta_prop_df, formulas)


if __name__ == "__main__":
    cli()
