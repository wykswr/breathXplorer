# breathXplorer

<!-- TOC -->
* [breathXplorer](#breathxplorer)
  * [Installation](#installation)
  * [Usage](#usage)
    * [File format](#file-format)
    * [Feature extraction](#feature-extraction)
    * [MS/MS spectra export](#msms-spectra-export)
    * [Feature alignment](#feature-alignment)
<!-- TOC -->

BreathXplorer is a swiss army knife for breath analysis. It provides a set of tools for breath analysis, including
feature extraction, feature alignment, and breath recognition.

## Installation

The package can be installed using pip, make sure you have python 3.7 to 3.9 installed.

`pip install breathXplorer`

## Usage

### File format
The input file should be mzML format. For the output, the single or aligned feature table can be exported as csv file.
The related MS/MS spectra can be exported as mgf file.

### Feature extraction
The feature extraction is performed using the `find_feature` function. The function takes as input the path to the mzMl file,
and returns an object containing the extracted feature table.

```python
from breathXplorer import find_feature

fs = find_feature("sample.mzML", False, .2, "Topological", 6)
```

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

### MS/MS spectra export
If you're using tandem MS, you can also export the MS/MS spectra as mgf file using the `to_mgf` function:

```python
from breathXplorer import retrieve_tandem

tandem = retrieve_tandem("sample.mzML", fs, 0.005)
tandem.to_mgf("ms2.mgf")
```

### Feature alignment
The `merge_feature` function takes as input a list of FeatureSet objects, and returns a Sample object.

```python
from breathXplorer import merge_result, find_feature

fss = [find_feature(f, False, .2, "Topological", 6) for f in ["sample1.mzML", "sample2.mzML", "sample3.mzML"]]
sample = merge_result(fss, ["sample1", "sample2", "sample3"])
```

Like the FeatureSet object, the Sample object contains the following information:

```python
sample.mz  # m/z values of the extracted features
sample.sample_name # the name of each sample
sample[96.7654] # get the total intensity with m/z value 96.7654 of all samples
sample.to_csv("aligned_table.csv") # export the feature table of all samples as csv file
```