[![DOI](https://zenodo.org/badge/360837483.svg)](https://zenodo.org/badge/latestdoi/360837483)

# Forgotten Books

## Paper
This repository holds the code and data accompanying the following paper:

"Forgotten Books: The Application of Unseen Species Models to the Survival of Culture" [2021], by Mike Kestemont, Folgert Karsdorp, Elisabeth de Bruijn, Matthew Driscoll, Katarzyna A. Kapitan, Pádraig Ó Macháin, Daniel Sawyer, Remco Sleiderink & Anne Chao.

> *Abstract*:  The study of ancient cultures is hindered by the incomplete survival of material artefacts, so that we commonly under-estimate the diversity of the cultural production in historic societies. To correct for this survivorship bias, we apply unseen species models from ecology and gauge the loss of narratives from medieval Europe, such as the romances about King Arthur. The obtained estimates are compatible with the scant historic evidence. Besides events like library fires, we identify the original evenness of populations as an overlooked factor in their stability in the face of immaterial loss. We link the elevated evenness in island literatures to parallel accounts of ecological and cultural diversity in insular communities. Our analyses call for a wider application of these methods across the heritage sciences.

## Code
The Jupyter notebooks under the `notebooks` folder hold all the Python code which we used for the analysis, including the additional experiments reported in the SI and the R code for reproducing our findings:
  - `analysis.ipynb`: code for the unseen species models
  - `analyze_pop.ipynb` (and `sim_pop.py`) for the evenness simulations
  - `geolocate.ipynb`: code used for plotting the heatmaps
  - `Copia.R`: the analogous R code to reproduce our findings

The code heavily relies on the open-source `copia` [package](https://github.com/mikekestemont/copia), co-developed by Mike Kestemont and Folgert Karsdorp, that is available from PyPI:

```bash
>>> pip install copia
```

The copia package is documented [here](https://copia.readthedocs.io/en/latest/). The code was executed in an Anaconda environment using Python 3.7.10. All dependencies can be installed from the `requirements.txt` file from the top-level directory in the repository:

```bash
>>> pip install -r requirements.txt
```


## Data
- The `datasets/master` folder contains spreadsheets (`.xlsx`) with the main work-document counts for the six main medieval vernaculars considered in the paper: Dutch, English, French, German, Icelandic, Irish (as well as Anglo-Norman). The compilation of these and the data format is detailed in the SI.
- The `datasets/geolocated` contains the data that was for plotting the heatmaps of the document dispersal (four vernaculars), with a latitude-longitude pair for each documents. Here, only documents are included that we were able to geolocalize approximately. If a document has multiple signatures (because its remnants are curated across multiple repositories), only the first signature is geolocalized.)


## Persistent archiving
Releases of this repository are sustainably mirrored on [Zenodo](https://zenodo.org/), ensuring long-term archival access to this material. Please consider citing the accompanying paper if you re-use this code for academic purposes.

## License
[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
