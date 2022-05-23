## Usage
(1) Get "sscount" file from mFold server
(2) Use `analyse.ipynb` (require Jupyter Notebook) to get statistical information
(3) Run `python open_close_scan.py` to get required ASO windows

- both (2) and (3) require modifying variable `sscount_file` to path of the target sscount file

## Requirement

```
python>=3.7
numpy
matplotlib
sklearn
xgboost
pandas
esmre
```