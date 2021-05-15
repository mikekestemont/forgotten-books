# Forgotten Books

## Paper
This repository holds the code and data accompanying the following paper:

"Forgotten Books: The Application of Unseen Species Models to the Survival of Medieval Literature" [2021], by Mike Kestemont, Folgert Karsdorp, Elisabeth de Bruijn, Matthew J. Driscoll, Katarzyna A. Kapitan, Pádraig Ó Macháin, Daniel Sawyer, Remco Sleiderink & Anne Chao.

> *Abstract*: Unseen species models are an important instrument in the unbiased monitoring of biodiversity. The applicability of these methods extends beyond ecology, however, as recent studies highlight ties between ecological and cultural (e.g., linguistic) diversity. Analogous to observation bias in biostatistics, the study of human culture is generally hindered by the incomplete survival of material artefacts. To help remedy issues of survivorship bias, we apply established estimators from ecology to statistically gauge the centuries-long loss of cultural diversity. We focus on the cultural domain of medieval narratives (heroic and chivalric literature), a hallmark of Western European culture, as attested in contemporary, handwritten documents. Our results confirm the overall severity of the sustained losses but reveal a wide range of survival rates across six representative vernaculars. Two insular literatures (Icelandic and Irish) combine a remarkably high evenness and survival rate, which we tentatively link to parallel observations about ecological diversity in island regions.

## Code
The Jupyter notebooks under the `notebooks` folder hold all the Python code which we used for the analysis, including the additional experiments reported in the SI:
  - `analysis.ipynb`: code for the analysis
  - `geolocate.ipynb`: code used for plotting the heatmaps

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
