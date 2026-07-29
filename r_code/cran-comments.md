# Resubmission

This is a resubmission of the `dynamicmultiplex` 1.0.0 submission, addressing
the review comments of 2026-07-29 (Konstanze Lauseker), together with the
1.1.0 feature release. Changes in response to the review:

* `\dontrun{}` removed. The single wrapped example
  (`animate_multilayer_gif()`) is now `\donttest{}`: it renders a GIF via
  'gganimate', which exceeds the 5-second example limit. All other examples
  are unwrapped and executable.

* No function writes to the user's filespace by default.
  `animate_multilayer_gif()` previously defaulted to writing
  "multilayer_animation.gif" in the working directory; `output_file` is
  now a required argument with no default, and the function errors with
  instructions (pointing to `tempfile()`) when it is missing. Its example
  writes only to `tempfile()`. No other function, example, or test writes
  outside `tempdir()`.

The version is 1.1.0 (rather than 1.0.1) because the resubmission also
includes the planned 1.1.0 release: an uncertainty-quantification overhaul,
cross-layer meta-community tracking, a Hungarian snapshot-matching method,
seed arguments for reproducible detection, and node-universe validation. See
NEWS.md.

## Test environments

* local Windows 11
* win-builder (devel and release)
* Ubuntu (R 4.3.3)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission (resubmission after review).

The possibly-misspelled words flagged in DESCRIPTION (Jaccard, Louvain,
interlayer, multislice) are standard network-analysis terms and method names.

## Downstream dependencies

There are currently no downstream dependencies for this package.
