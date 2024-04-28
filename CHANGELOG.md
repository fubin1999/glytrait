# Changelog

Currently, GlyTrait is in development and no stable version has been released yet. 
The following is a list of changes that have been made to the project.

## v0.1.4

**Added:**

- Sialic acid linkage configuration in the Streamlit app.
- A spinner is shown now when the Streamlit app is running.

**Fixed:**

- Fix the bug that the "mode" configuration was ignored by the Streamlit app.

## v0.1.3

**Added:**

- The Streamlit app!
- A new method for the `Experiment` class: `Experiment.run_workflow()`.

## v0.1.2

**Added:**

- `Experiment` class is added as the new api of GlyTrait, instead of the old `GlyTrait`. 
  This might be part of the python package api in the future.

**Changed:**

- The formula file format is much simpler: any non-blank lines not starting with
  "#" will be parsed as a formula. 
  The formula descriptions are converted into comments (starting with "#").
- `GlyTraitInputData` now checks the validity of data.
  This allows dynamic change the input data after a `GlyTraitInputData` instance
  is created.

**Removed:**

- `GlyTrait` class is removed.

## v0.1.1

**Changed:**

- When provided a custom formula file, the built-in formula file is no longer used.

## v0.1.0

Initial release of GlyTrait.
