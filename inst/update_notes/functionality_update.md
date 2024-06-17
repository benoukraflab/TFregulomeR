## Current Functionalities v2.3.3 - 17 June 2024

Along with an update to the function _greatAnnotate_ we have added the ability to run the GREAT analysis locally.

## Previous  Version
### Functionalities version 2.3.2 - 3 June 2024
1) Retired rbokeh HTML output due to the package no longer being maintain.
2) Recreated the HTML output of _greatAnnotate_ using the three new packages: plotly, crosstalk, and DT.

### Functionalities version 2.3.1 - 26 April 2024
Removed deprecated format for _aes_ data in ggplot2.

### Functionalities version 2.3.0 - 16 April 2024
bug fixed allowing for backwards compatibility with older local database versions.

### Functionalities version 2.2.3 - 8 March 2024
bug fixed within _.searchMotif_.

### Functionalities version 2.2.2 - 29 February 2024
1) Fixed export error for _apiRequest_.
2) Exported _formMatrixFromSeq_ for the [forkedTF](https://github.com/benoukraflab/forkedTF) R library.

### Functionalities version 2.2.1 - 23 January 2024
1) Built a class to handle all server side or local calls to the Tregulome compendium.
2) Replaced all deprecated compendium calls with a class instantiation and/or method call to the new _apiClass_.

### Functionalities version 2.2.0 - 19 December 2023
1) Overhauled the MEME and TRANSFAC functionality throughout the package.
2) Modularized _intersectPeakMatrix_ to allow for integration within the [forkedTF](https://github.com/benoukraflab/forkedTF) R library.
3) Exported previously hidden helper functions for the [forkedTF](https://github.com/benoukraflab/forkedTF) R library.

### Functionalities version 2.1.2 - 24 May 2023
bug fixed within _intersectPeakMatrix_ ability to save MethMotif logos.

### Functionalities version 2.1.1 - 4 May 2023
Bug fixed in _local\_db\_path_ argument.

### Functionalities version 2.1.0 - 24 April 2023
1) Added new feature to call a local version of the TFregulome compendium allowing for faster execution speeds.
2) Functionality rolled out as argument _local\_db\_path_ to _commonPeaks_, _dataBrowser_, _exclusivePeaks_, _genomeAnnotate_, _greatAnnotate_, _intersectPeakMatrix_, _loadPeaks_, _motifDistrib_, _searchMotif_, and _toTFBSTools_.
3) The local version of TFregulome's database can be found [here](https://methmotif.org/API_TFregulomeR/downloads/)

### Functionalities version 2.0.2 - 17 April 2022
Changed to Canada server as default (Singapore server retired)

### Functionalities version 2.0.1 - 7 May 2021
bug fixed regarding http failure

### Functionalities version 2.0.0 - 5 March 2020
1) Added new TFregulomeR compendium server in Canada. Users now can choose servers between Singapore and Canada.
2) linked to mouse TFregulomeR compendium in Canada server, coming soon in Singapore server.

### Functionalities version 1.2.2 - 10 January 2020
1) Modified the functions _intersectPeakMatrix_ and _intersectPeakMatrixResult_ to allow profiling of users' input external signals during pair-wise comparison analysis.
2) Added a new function _interactome3D_ which could directly generate a dynamic 3D interface report showing TF interactome coupled with CpG methylation and users' input external signal values.

### Functionalities version 1.2.1 - 25 October 2019
1) Improved the function _greatAnnotate_, which now allows to use peak regions in hg38 directly with rGREAT >= 1.16.1 (no need liftover package).

### Functionalities version 1.2.0 - 6 September 2019
1) Changed the function name _TFBSBrowser_ to _dataBrowser_.
2) More detailed information about motifs, such as consistency with exisiting database and correlation between different motif callers, added in the output of _dataBrowser_.
3) Fixed some bugs in function _cofactorReport_.
4) Added y axis and number of peaks with motif in the output of _plotLogo_.

### Functionalities version 1.1.0 - 31 July 2019

1) Browse the TFregulomeR data warehouse (TFBSBrowser)
2) Load TF peaks (loadPeaks)
3) Search motif matrix and DNA methylation score matrix (searchMotif)
4) Plot motif or MethMotif logo (plotLogo)
5) Export motif matrix and DNA methylation score matrix (exportMMPFM)
6) Get context-independent peaks along with DNA methylation profiles (commonPeaks & commonPeakResult)
7) Get context-dependent peaks along with DNA methylation profiles (exclusivePeaks & exclusivePeakResult)
8) Form a intersected matrix between two lists of peak sets along with DNA methylation profiles and read enrichments, for interactome and co-binding partner studies (intersectPeakMatrix & intersectPeakMatrixResult)
9) Automatically generate a PDF report for TF co-factors along with motif sequences, DNA methylation (within motif and in 200bp regions) and read enrichments (cofactorReport).
10) Plot the TFBS distribution in a given list of peak sets (motifDistrib & plotDistrib)
11) Annotate peak genomic locations (genomeAnnotate)
12) Annotate ontologies of target genes by a peak set (greatAnnotate)
13) Convert a motif matrix to a PFMatrix calss object for *TFBSTools* package (toTFBSTools)

### Functionalities version 1.0.0 - 7 May 2019
1) Browse the TFregulomeR data warehouse (TFBSBrowser)
2) Load TF peaks (loadPeaks)
3) Search motif matrix and DNA methylation score matrix (searchMotif)
4) Plot motif or MethMotif logo (plotLogo)
5) Export motif matrix and DNA methylation score matrix (exportMMPFM)
6) Get context-independent peaks along with DNA methylation profiles (commonPeaks & commonPeakResult)
7) Get context-dependent peaks along with DNA methylation profiles (exclusivePeaks & exclusivePeakResult)
8) Form a intersected matrix between two lists of peak sets along with DNA methylation profiles, for interactome and co-binding partner studies (intersectPeakMatrix & intersectPeakMatrixResult)
9) Plot the TFBS distribution in a given list of peak sets (motifDistrib & plotDistrib)
10) Annotate peak genomic locations (genomeAnnotate)
11) Annotate ontologies of target genes by a peak set (greatAnnotate)
12) Convert a motif matrix to a PFMatrix calss object for *TFBSTools* package (toTFBSTools)
