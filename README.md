# shortRes
A Sage package for computing the short resolution of a weighted homogeneous ideal. For more details about the functions included, please check the associated paper in [arxiv doi] 

Examples can be found in the help section of each function (function? will return this section). 

**How to use**: download shortRes.py and write `load("shortRes.py")` in Sage. Depending on your current directory, this may require specifying the path of the package, e.g., `load("/home/user/shortRes.py")`.  

**Small example**:
```python
load("GHWs.py")
C = codes.BinaryReedMullerCode(1,4)
hierarchy(C)
```
Output: [8, 12, 14, 15, 16]


**List of functions**: We provide several main functions, as well as some auxiliary functions for working with generalized Hamming weights. 
  - **Main functions**: 'GHW', 'hierarchy', 'RGHW', 'rhierarchy', 'higher_spectrum', 'rhigher_spectrum', 'matrix_supp', 'subspaces', 'num_subspaces', 'wei_duality'.
  - **Low memory main functions**: 'GHW_low_mem', 'hierarchy_low_mem', 'RGHW_low_mem', 'rhierarchy_low_mem', 'higher_spectrum_low_mem', 'rhigher_spectrum_low_mem'.
  - **Auxiliary functions**: 'vecwt', 'coltw', 'standard', 'is_cyclic', 'bch_bound', 'information'.

**Tests**: it is possible to test that the functions are working propertly by running the test test_GHWs.sage (test_GHWs_low_mem.sage for the low memory functions). This can be done by writing `sage test_GHWs.sage`. This assumes the folder structure follows that of this repository. Otherwise, the line `load(path2)` from the test file has to be changed to specify the path of GHWs.py. The test should take between 130s and 500s, depending on whether the low memory functions are used or not and the processor's performance. The rest of the files are performance tests used to obtain the tables and graphs of the associated paper.

**Citation**: if you use this implementation for your research, please consider citing the paper https://doi.org/... and/or this repository.
