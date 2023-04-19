# breathXplorer

BreathXplorer is a toolkit for the analysis of breath metabolomics data. 

## Installation

`pip install breathXplorer`

## Usage

### File format
The input file should be mzML format.

### Feature extraction
The feature extraction is performed using the `find_feature` function. The function takes as input the path to the mzMl file,
and returns a pandas dataframe with the extracted features.

```python
from breathXplorer import find_feature

tb = find_feature("sample.mzML", False, .2, "Topological", 6)
```

The extracted feature table has the m/z values as index, and the total intensity over time as the first column,
the rest of the columns are times.

### Merge feature tables
The `merge_feature` function takes as input a list of extracted feature table, and returns a merged feature table.

```python
from breathXplorer import merge_result, find_feature

tbs = [find_feature(f, False, .2, "Topological", 6) for f in ["sample1.mzML", "sample2.mzML", "sample3.mzML"]]
tb = merge_result(tbs, ["sample1", "sample2", "sample3"])
```

The merged feature table has the m/z values as index, and sample names as columns. The values are the total intensity over time.