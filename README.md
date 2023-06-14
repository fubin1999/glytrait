![logo](img/logo.png)

# GlyTrait

## Overview

Glycan derived trait is a more insightful way to analysis glycomics data. However, currently there
lacks a tool for automatically calculating derived traits, while mannual calculating is
cumbersome, time-consuming and error-prone. GlyTrait is a tool designed for calculation N-glycan
traits merely from abundance information and glycan structures.

## Installation

### Requirement

```
python >= 3.10
```

If python hasn't been install, download it from
its [its website](https://www.bing.com/search?q=python&form=APMCS1&PC=APMC), or
use [Anaconda](https://www.anaconda.com/download/) if you like.

### Using pipx (recommended)

pipx is a tool to help you install and run end-user applications written in Python. It's roughly
similar to macOS's brew, JavaScript's npx, and Linux's apt.

#### Insatll pipx

Install pipx following its [Document](https://pypa.github.io/pipx/installation/).

#### Install GlyTrait from PyPi

```shell
pipx install glytrait
```

#### Or install GlyTrait locally

1. Clone the repository (or download it manually):

```shell
git clone https://github.com/fubin1999/glytrait.git
```

2. Move to the repository:

```shell
cd glytrait
```

3. Install glyTrait with pipx:

```shell
pipx install .
```

## Usage

### Quick start

```shell
glytrait data.csv
```

That's it! If everything goes well, a xlsx file named "data_glytrait.xlsx" will be saved to the
same directory with data.csv. Inside the xlsx file are three sheets, the trait values for each
sample, the trait descriptions, and the meta properties GlyTrait generated to calculated derived
traits. You don't need the third one for common situations, but a peek of it might help you better
understand how GlyTrait works.

### Input file format

The input csv file should have at least 3 columns:

1. The **first** column with header "Glycan ID" is the identifier of the glycans. Glycan ID could
   theoretically be any string, but we recommend using condense compostion for clarity, e.g. "
   H5N2", "H5N4F1S2".
2. The **second** column with header "Structure" is the structure string of the glycans. We only
   support condensed GlycoCT now as it can be exported from GlycoWorkbench. A complete structure
   with all linkage specified is not necessary, because GlyTrait actually use the topology
   properties of glycan structures regardless of the linkage information. However, sialic acid
   linkage is needed if you want to calculate the sialic-acid-linkage traits. See below for more
   information.
3. **From the third column on** is the glycan abundance for different samples, with sample names
   as the headers. Normalization is not necessary because GllyTrait will carry out a Total
   Abundance Normalization for all samples before calculating derived traits.

**Note that Glycan ID and Structures should all be unique!**

An example input file would be like:

| Glycan ID | Structrue | Sample1 | Sample2 | Sample3 |
|-----------|-----------|---------|---------|---------|
| H3N3F1    | RES...    | 0.0417  | 0.0503  | 0.0354  |
| H3N4      | RES...    | 0.0233  | 0.0533  | 0.0593  |
| H3N4F1    | RES...    | 0.0123  | 0.0133  | 0.0194  |

This file contains 3 glycans (H3N3F1, N3N4 and N3N4F1) and three samples (Sample 1, 2, and 3).

### Sialic-acid-linkage traits

Sialic acids can have different linkages for N-glycans (e.g. α2,3 and α2,6). Different sialic acid
linkage might have different biological functions. GlyTrait supports calculating derived traits
regarding to these linkages. To get these traits, you have to provide sialic acid linkage
information in the Structures for **all** glycans in the input file, and use the `-s`
or `--sia_linkage` flag:

```shell
glytrait data.csv -s
```

### Specified output path

You might noticed before that GlyTrait save the output file to the same directory as the input
file with a "_glytrait" suffix. You can specify the output file path by using the "-o" or "
--output_file" option:

```shell
glytrait data.csv -o output.xlsx
```

Note that a ".xlsx" suffix is needed.

### Filtering

Some derived traits might not be useful for your analysis. For example, some traits might all have
the same value for all samples, or be NaN due to zero being in the denominator. GlyTrait rules out
these traits by default. If you want to keep these traits, use the "--no-filter" option:

```shell
glytrait data.csv --no-filter
```

### Hypothesis testing

GlyTrait supports hypothesis testing for direct and derived traits. To use this feature, you need
to provide a csv file containing the sample grouping information. The csv file should have two
columns, the first column is the sample name, and the second column is the group name. An example
group file would be like:

| Sample  | Group  |
|---------|--------|
| Sample1 | Group1 |
| Sample2 | Group1 |
| Sample3 | Group2 |
| Sample4 | Group2 |

Then use the "-g" or "--group_file" option to specify the group file:

```shell
glytrait data.csv -g group.csv
```

GlyTrait will carry out Mann-Whitney U Test for two groups, and Kruskal-Wallis H Test for more
than two groups. The Benjamini-Hochberg procedure is used to correct the p-values. For 
more-than-two-groups situations, the Mann-Whitney U Test will be used for post-hoc test. GlyTraits
uses non-parametric tests for the sake of robustness.

### The GlyTrait Formula (Advanced)

Currently 164 derived traits (including sialic acid linkage traits) is included in the GlyTrait
tool. This list is curated from literature and covers nearly all derived traits reported. However,
you may want to add new traits to fulfill your own need. GlyTriat using a meta-properties-oriented
formula system, the **GlyTrait Formula**, to represent the meaning of a trait. **GlyTrait Formula
** is a versatile tool allowing you to add your own traits.

#### GlyTrait Formula overview

**GlyTrait Formula** is a text based formula representing sys allowing you to combine various meta
properties to get new traits. As an overview, the formula should be in the format of:

![formula](img/formula.png)

- <Name> is the name of the formula, which will be used as the column name in the output file.
- <Numerator1>, <Numerator2>, ... are the meta properties of the numerator.
- <Denominator1>, <Denominator2>, ... are the meta properties of the denominator.
- () is necessary to group the numerator and denominator for clearness, even if there is only one
  meta property on the numerator and the denominator.
- All numerators and denominators should be separated by "*". If there is only one numerator or
  denominator, "*" is not needed.
- The "/" is necessary to separate the numerator and denominator. If the trait is a proportion of
  something within something, use "//" instead of "/".
- Some formulas may use all glycans as the denominator. In this case, use (.) as the denominator.
  A combination of "." and other meta properties (e.g. (. * isComplex)) is not
  allowed.
- <Coefficient> is the coefficient multiplied to the final result. The coefficient part (* <Coefficient>) is optional. If omitted, the coefficient is assumed to be 1.

For more details and step-by-step instructions on how to write a formula, refer to the template
file generated by `glytrait -t some_dir` (see setion below).

#### Get the formula template file

To start using **GlyTrait Formula**, you first need to get a formula template file by using the "
-t" or "--save_template" option:

```shell
glytrait -t some_dir
```

An template txt file will be saved in the given "some_dir" directory. "some_dir" needed to exist,
GlyTrait will create this folder for you if it doesn't exist. A detailed instruction on how to
write trait formulas is in that file.

#### Using the custom formulas

Once you have editing the template file with new formulas in, you can try them out. For
incoperating the custom formulas, use the "-f" or "--formula_file" option:

```shell
glytrait -f custom_formulas.txt
```

That will do, GlyTrait will calculate the default traits, plus your cumtom traits.
A few things to be noted:

- If a custom formula has the same name with a default formula, it will be ignored.
- If more than one custom formulas have the same name, only the first one will be used.
- If any formula is in wrong format, GlyTrait will raise an error to inform you.

## License

[MIT License](LICENSE)
