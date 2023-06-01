# breathXplorer

[![PyPI](https://img.shields.io/pypi/pyversions/breathXplorer)](https://pypi.org/project/breathXplorer/)
[![Python package](https://github.com/wykswr/breathXplorer/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/wykswr/breathXplorer/actions/workflows/python-package.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- TOC -->

* [breathXplorer](#breathxplorer)
    * [Introduction](#introduction)
    * [Quick start](#quick-start)
        * [Installation](#installation)
        * [File format](#file-format)
        * [Feature extraction](#feature-extraction)
        * [Feature alignment](#feature-alignment)
    * [Utilities](#utilities)
        * [MS/MS spectra export](#msms-spectra-export)
        * [Peak recognition](#peak-recognition)
    * [Citation](#citation)

<!-- TOC -->

## Introduction

BreathXplorer is a swiss army knife for breath analysis. It provides a set of tools for breath analysis, including
feature extraction, feature alignment, and breath recognition.

## Quick start

### Installation

The package can be installed using pip, supported python version are 3.7, 3.8, 3.9 and 3.10.
It's recommended to install the package in a virtual environment or a conda environment,
after activating the environment, run the following command to install the package:

`pip install breathXplorer`

### File format

The input file should be in mzML format. For the output, the single or aligned feature table can be exported as csv
file.
The related MS/MS spectra can be exported in an mgf file.

### Feature extraction

Feature extraction is used to find the volatile organic compound (VOC) in the breath sample.
The feature extraction is performed using the `find_feature` function. The function takes the path to an mzMl file
as input, and returns an object containing the extracted feature table:

```python
from breathXplorer import find_feature

fs = find_feature("sample.mzML", False, .8, "Topological", 6)
```

The `False` indicates that the input file is not a line spectrum, and the `.8` controls the quality of the extracted
features, higher value means higher quality. The `"Topological"` indicates the algorithm used for feature extraction,
the other option is `"Gaussian"`. The `6` is the prior knowledge of the number of breath in an experiment (only used
for Gaussian algorithm).

The `fs` is a FeatureSet object, it contains the following information:

```python
fs.mz  # m/z values of the extracted features
fs.scan_time  # scan time of the experiment
fs.intensity  # the total intensity of each feature (calculated by integrating the intensity over scan time)
len(fs)  # the number of extracted features
fs[96.7654]  # get the intensity with m/z value 96.7654 over scan time
```

FeatureSet object can be exported as csv file using the `to_csv` method:

```python
fs.to_csv("feature_table.csv")
```

### Feature alignment

The `merge_feature` function takes as input a list of FeatureSet objects, and returns a Sample object. It aligns the
features with the similar m/z value from different samples, and calculate the total intensity of each feature in each
sample. To use the function, you can do the following:

```python
from breathXplorer import merge_result, find_feature

fss = [find_feature(f, False, .8, "Gaussian", 6) for f in ["sample1.mzML", "sample2.mzML", "sample3.mzML"]]
sample = merge_result(fss, ["sample1", "sample2", "sample3"])
```

The first statement creates a list of FeatureSet objects, the second statement aligns those FeatureSet objects.
The `fss` is a list of FeatureSet objects, the `["sample1", "sample2", "sample3"]` is the customizable names assigned
to each sample.

Just like the FeatureSet object, the Sample object contains the following information:

```python
sample.mz  # m/z values of the extracted features
sample.sample_name  # the name of each sample
len(sample)  # the number of extracted features
sample[96.7654]  # get the total intensity with m/z value 96.7654 of all samples
sample.to_csv("aligned_table.csv")  # export the feature table of all samples as csv file
```

## Utilities

### MS/MS spectra export

If you're using tandem MS, you can also export the MS/MS spectra as mgf file using the `to_mgf` function:

```python
from breathXplorer import retrieve_tandem

tandem = retrieve_tandem("sample.mzML", fs, 0.005)
tandem.to_mgf("ms2.mgf")
```

The file `sample.mzML` contains tandem MS data. The `fs` object represents either a FeatureSet or a Sample, from which
the
m/z values are extracted. These m/z values are then used to obtain the corresponding MS/MS spectra. A tolerance
of `0.005`
is applied to the m/z values. Any MS/MS spectra with a difference smaller than `0.005` between the precursor m/z value
and
the feature's m/z value will be retrieved.

### Peak recognition

The total ion current (TIC) of mass spectral data from breath analysis exhibits two key properties that preclude the use
of existing peak detection methods but facilitate the development of specialized breath analysis algorithms:

1. Breath peaks are characteristically irregular and contain a substantial number of subsidiary peaks.
2. Prior knowledge about the expected number of breath peaks is typically available.

We developed 2 algorithms to detect breath peaks from the TIC of mass spectral data from breath analysis: Topological
and
Gaussian mixture model (GMM). Both the 2 methods share the similar interface, and can be imported as follows:

```python
from breathXplorer.find_peak import find_peak, find_gaussian
```

Both function takes in two arrays, `x` representing the time points and `y` representing the corresponding values such as
concentration or intensity. Besides, when using `find_gaussian`, the parameter `n` determines the number of peaks to be considered. The output is a list of
tuples, where each tuple contains the start and end time point values of the identified ranges, along with the maximum
intensity values of the respective peaks.

## Citation

If you use this package in your research, please cite the following paper:

