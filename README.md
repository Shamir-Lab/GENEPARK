# To get the results of the GENEPARK project:
1. Download all code and data files into the same directory
2. Open reproduce_geneparks_results.R 
3. Set the working dir to the directory in which all files are in
4. Run the analyses to get the models, signature, and figures. Note that some analyses are slow (especially the fSVA-based).

# The gene expression data files 
These are all the compressed .7z files. Here are important comments on these files:
1. The training set files should be merged  - they were split because of the file size (>25mb). These are the training set samples after preprocessing.
2. The training and validation combined data files should be merged  - they were split because of the file size (>25mb). These are the training and validation combined set after preprocessing. These data were used for the training set and for obtaining the biomarker.
3. Validation set - a set of samples that were used for initial tests and for tuning.
4. The final test set

License
=======

* Low level code for dealing with binary encoded k-mers and strings, and for Bloom filters is derived from the original minia implementation, http://minia.genouest.org/; these components, mostly unmodified, are distributed under a GPL 3.0 license

* Original code is distributed under the BSD 3 clause license.
