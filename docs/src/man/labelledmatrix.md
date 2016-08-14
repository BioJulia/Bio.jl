# Reference: Bio.LabelledMatrices - Labelled square matrix parsers

```@meta
CurrentModule = Bio.LabelledMatrices
```

This module contains a function `readlsm`, that reads labelled square matrices.

Such matrices are commonly used to store pairwise measures of similarity or
dissimilarity between two sets. This format consists of a delimited file,
with text labels on the upper and left sides of a square matrix. Such matrices
are often but not exclusively symmetric.

An example would be: 

```
    A   B
A   1   2
B   2   1
```

(Note that we use spaces as a delimiter, where tabs would normally be used)

