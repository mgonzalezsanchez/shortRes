# shortRes
A Sage package for computing the short resolution of a weighted homogeneous ideal. For more details about the functions included, please check the associated paper in [arxiv doi] 

Examples can be found in the help section of each function (function? will return this section). 

**How to use**: download shortRes.py and write `load("shortRes.py")` in Sage. Depending on your current directory, this may require specifying the path of the package, e.g., `load("/home/user/shortRes.py")`.  

**Small example**:
```python
load("shortRes.py")
A = numpy.array([[7,2,3],[1,8,3],[3,8,1],[12,0,0],[0,12,0],[0,0,12]]);
r,I = toricIdeal(A,QQ);
res = shortRes(r,I,1);
```
Output: 

 ![Output](betti.PNG)

**List of functions**: We provide several main functions, as well as an auxiliary function. 
  - **Main functions**: 'toricIdeal', 'shortRes', 'shortResDim3', 'pruningDim3'.
  - **Auxiliary function**: 'sortLex'.

**Citation**: if you use this implementation for your research, please consider citing the paper https://doi.org/... and/or this repository.
