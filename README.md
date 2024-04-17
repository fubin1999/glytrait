![logo](img/logo.png)

# GlyTrait

[![PyPI - Version](https://img.shields.io/pypi/v/glytrait)](https://pypi.org/project/glytrait/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/glytrait)
[![GitHub License](https://img.shields.io/github/license/fubin1999/glytrait)](https://github.com/fubin1999/glytrait/blob/main/LICENSE)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e3e6a19dccb749f786264247738fa585)](https://app.codacy.com/gh/fubin1999/glytrait/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

## Overview

*Q: What is GlyTrait?*

*A: GlyTrait is a tool to calculate derived traits for N-glycans.*

*Q: Sounds cool! So..., what are derived traits again?*

*A: Well, derived traits are artificial variables that summarize certain aspects of 
the glycome. 
For example, the number of antennas, the number of fucoses, the number of sialic acids,
and their interactions, etc.*

## Contents

- [Installation](#installation)
    - [Requirement](#requirement)
    - [Using pipx (recommended)](#using-pipx-recommended)
- [Usage](#usage)
    - [Quick start](#quick-start)
    - [Options](#options)
    - [Mode](#mode)
    - [Input file format](#input-file-format)
    - [Specify output path](#specify-the-output-path)
    - [Preprocessing](#preprocessing)
    - [Sialic-acid-linkage traits](#sialic-acid-linkage-traits)
    - [Post-filtering](#post-filtering)
- [License](#license)

## Installation

### Requirement

```
python >= 3.10
```

If python hasn't been installed, 
download it from [its website](https://www.bing.com/search?q=python&form=APMCS1&PC=APMC),
or use [Anaconda](https://www.anaconda.com/download/) if you like.

### Option 1: Using pipx (recommended)

`pipx` is a tool to help you install and run end-user applications written in Python.
It's roughly similar to macOS's brew, JavaScript's npx, and Linux's apt.

#### Install pipx

Install pipx following its [Document](https://pypa.github.io/pipx/installation/).

#### Install GlyTrait from PyPi

```shell
pipx install glytrait
```

### Option 2: Using pip

#### 1. Create a new conda environment (optional)

```shell
conda create -n glytrait python=3.10
```

and activate the environment.

```shell
conda activate glytrait
```

#### 2. Install GlyTrait from PyPi

```shell
pip install glytrait
```

## Usage

### Quick start

1. Download the [example files](https://github.com/fubin1999/glytrait/tree/main/example_data)
to a new directory: `glytrait_example`.

2. Open a terminal and navigate to the directory:

```shell
cd glytrait_example
```

3. Run the following command:

```shell
glytrait abundance.csv structures.csv
```

That's it! 
If everything goes well, a folder named "abundance_glytrait" will be created
in the same directory with the abundnce.csv file.
Inside the directory are four files:

1. `derived_traits.csv`: all derived traits calculated by GlyTrait.
2. `derived_traits_filtered.csv`: derived traits after post-filtering 
   (see [Post-filtering](#post-filtering)).
3. `glycan_abundance_processed.csv`: the glycan abundance after preprocessing
   (see [Preprocessing](#preprocessing)).
4. `meta_properties.csv`: the meta-properties of all glycans.

The detailed format of the input file will be introduced in the
[Input file format](#input-file-format) section.

### Options

*This section is intended to give an overview of the CLI interface.
Feel free to skip it right now.*

As a glance, GlyTrait supports the following options:

|          Option           | Description                                                                                                     |
|:-------------------------:|:----------------------------------------------------------------------------------------------------------------|
|          --help           | Show the help message and exit.                                                                                 |
|        -m, --mode         | The mode. "S" or "structure" for structure mode, "C" or "composition" for composition mode. Default: structure. |
|       -o, --output        | The output path. Default: the same directory with the input file.                                               |
| -r, --filter-glycan-ratio | The proportion of missing values for a glycan to be ruled out. Default: 0.5.                                    |
|    -i, --inpute-method    | The imputation method. "min", "mean", "median", "zero", or "lod". Default: "min".                               |
|   -c, --corr-threshold    | The correlation threshold for collinearity filtering. Default to 1.0.                                           |
|     -l, --sia-linkage     | Flag to include the sialic acid linkage traits.                                                                 |
|        --no-filter        | Flag to turn off post-filtering of derived traits.                                                              |
|        -g, --group        | The group file.                                                                                                 |
|    -f, --formula-file     | The custom formula file to use.                                                                                 |
|  -b, --builtin-formulas   | The directory path to save a copy of the built-in formulas.                                                     |

The following sections will introduce these options in detail.

### Mode

GlyTrait has two modes: the **"structure"** mode and the **"composition"** mode.
In the "structure" mode,
GlyTrait will calculate derived traits based on the topology properties of glycan structures.
In the "composition" mode, GlyTrait will make educated guesses on the structure properties 
based on the glycan composition.

Note that the "composition" mode has uncertainties to some extent. Specifically:

1. Estimating the number of Gal based on composition is not possible for hybrid glycans,
   so GlyTrait will calculate the number of Gal assuming there are no hybrid glycans.
   kily, hybrid glycans are usually in low abundance,
   so the algorithm is a good approximation for most cases.
2. Estimating the number of branches is not possible based on composition,
   so GlyTrait will roughly classify glycans into two categories: 
   low-branching and high-branching.
   Glycans with N > 4 (including bisecting diantenary glycans) are considered as high-branching,
   while others as low-branching.
3. Telling hybrid glycans from mono-antenary complex glycans is not possible based on composition,
   so GlyTrait will not classify glycans into complex, hybrid and high-mannose.

Due to the ambiguities above, we recommend using the "structure" mode if possible.

You can specify the mode by the "-m" or the "--mode" option:

```shell
glytrait abundance.csv composition.csv -m composition
```

Or in short:

```shell
glytrait abundance.csv composition.csv -m C
```

The default mode is the "structure" mode, as in the quick start example. 
Thus, using `glytrait -m structure` or `glytrait -m S` is equivalent to `glytrait` alone.
If you might use both modes in a project, 
we recommend using the "-m" option to specify the mode explicitly.

### Input file format

At least two files are needed for GlyTrait to work:

#### 1. The abundance file

A csv file with samples as rows and glycan IDs as columns.
An example file would be like:

| Sample  | Glycan1 | Glycan2 | Glycan3 |
|---------|---------|---------|---------|
| Sample1 | 0.0417  | 0.0503  | 0.0354  |
| Sample2 | 0.0233  | 0.0533  | 0.0593  |
| Sample3 | 0.0123  | 0.0133  | 0.0194  |

The header of the first column should be "Sample",
and the header of the other columns should be glycan IDs.
Glycan IDs can be any string, e.g., the composition strings ("H3N4").

**Both glycan IDs and samples should be unique.**

#### 2. The structure file (or the composition file)

A csv file with two columns: "GlycanID" and "Structure" (or "Composition").
An example file would be like:

| GlycanID | Structure |
|----------|-----------|
| Glycan1  | RES...    |
| Glycan2  | RES...    |
| Glycan3  | RES...    |

The "GlycanID" column should contain all glycan IDs in the abundance file.
The "Structure" column should contain the structure strings of the glycans.
For now, only the GlycoCT format is supported.
In the "composition" mode, the second column should be "Composition" instead of "Structure",
and the composition strings should be used instead of the structure strings.
Condensed format ("H3N4F1S1") is supported.

### Specify the output path

You might have noticed before that GlyTrait saves the output file to the same directory 
as the abundance file with a "_glytrait" suffix.
You can specify the output file path by using the "-o" or "--output-file" option:

```shell
glytrait abundance.csv structure.csv -o output
```

### Preprocessing

GlyTrait will carry out a preprocessing step before calculating derived traits.
The following steps will be done:

- Remove glycans with missing values in more than a certain proportion of samples.
- Impute missing values.
- Perform Total Abundance Normalization.

In the glycan-filtering step, the proportion threshold could be specified by the "-r" or the
"--filter-glycan-ratio" option.
The default value is 1, which means no glycan will be removed.
You can change this value to 0.5 by:

```shell
glytrait abundance.csv structure.csv -r 0.5
```

The imputation method could be specified by the "-i" or the "--impute-method" option.
The default method is "zero",
which means missing values will be imputed by 0.
Other supported methods are "mean", "median", "zero", "lod".
You can change the imputation method to "min" by:

```shell
glytrait abundance.csv structure.csv -i min
```

A full list of supported imputation methods are:

- "min": impute missing values by the minimum value of a glycan within all samples.
- "mean": impute missing values by the mean value of a glycan within all samples.
- "median": impute missing values by the median value of a glycan within all samples.
- "zero": impute missing values by 0.
- "lod": impute missing values by the limit of detection (LOD) of the equipment. The LOD of a
  glycan is defined as the minimum value of the glycan within all samples divided by 5.

### Sialic-acid-linkage traits

Sialic acids can have different linkages for N-glycans (e.g., α2,3 and α2,6).
Different sialic acid linkage has different biological functions.
GlyTrait supports calculating derived traits regarding these linkages.
To use this feature, you need to have siaic acid linkage information.

In the structure mode, the "Structure" column or the structure file should contain the linkage
information.
Only linkage information about sialic acids is needed.
This can be easily done using GlycoWorkbench.

In the composition mode, the "Composition" column must contain the linkage information.
GlyTrait uses a common notation for sialic acid with different linkages:
"E" for a2,6-linked sialic acids, and "L" a2,3-linked sialic acids.
For example, "H5N4F1E1L1" contains two sialic acids, one is a2,6-linked, 
and the other is a2,3-linked.

You can use the "-l" or "--sia-linkage" option to include sialic-acid-linkage traits:

```shell
glytrait abundance.csv structure.csv -l
```

Note that if you use this option, all glycans with sialic acids should have linkage information.
That is to say, all structure strings should have structure information in the structure mode and
no "S" in composition strings in the composition mode.

### Post-Filtering

Not all derived traits are informative.
For example, some traits might have the same value for all samples.
Some traits might be highly correlated with others.

GlyTrait carries out a two-step post-filtering process to remove these uninformative traits.
First, traits with the same value for all samples will be removed.
Second, highly correlated traits will be pruned,
keeping only the traits considering more glycans.

GlyTrait filters out highly correlated traits, using a "trait family tree" filtering method.
Briefly, for two correlated traits, the "parent" trait,
which normally considers more glycans, will be kept.
For example, for the two high-correlated traits: A2FG and A2G, the latter will be kept 
because it is more general, and more robust for considering more glycans.
Thanks to the dynamic "trait family tree" generated by GlyTrait,
user-defined traits will also be considered in this filtering process.

By default, GlyTrait only filters trarits with Pearson correlation coefficient of 1,
i.e., traits with perfect collinearity.
This threshold can be changed by the "-c" or "--corr-threshold" option:

```shell
glytrait abundance.csv structure.csv -c 0.9
```

Setting the threshold to -1 will turn off the colinearity filtering. 
To turn off postfiltering all together, use the "--no-filtering" option:

```shell
glytrait abundance.csv structure.csv --no-filtering
```

## License

[MIT License](LICENSE)
