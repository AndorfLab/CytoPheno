
# CytoPheno 

## Overview
This is a R Shiny application that intakes post-clustered cytometry data and denotes marker description patterns and descriptive cell type names.

Specifically, the full application takes post-clustering flow or mass cytometry expression data as the input, returns marker definitions per cluster (Part 1), standardizes marker names (Part 2), and finally matches to cell type names (Part 3). 

![Figure 1](https://github.com/user-attachments/assets/fd6ec134-f358-43cd-b891-cce9790c5bf2)


Post-clustered expression data is inputted via the **Input Expression Data** tab within the application. If the user prefers to directly input marker descriptions (e.g. CD4+, CD8-) and skip Part 1, they can do so through the **Input Marker Descriptors** tab.


## Index
**1. [Installation](https://github.com/AndorfLab/CytoPheno/wiki/1.-Installation)** – Installation instructions and R packages that are required to run the application  

**2. [Example Data](https://github.com/AndorfLab/CytoPheno/wiki/2.-Example-Data)** – Example datasets that can be used to demonstrate the application   

**3. Input Expression Data Tab** – Use when uploading post-clustered expression data into the application (Parts 1-3)  

* **3.1. [Input Expression Data and Default References](https://github.com/AndorfLab/CytoPheno/wiki/3.1.-Input-Expression-Data-and-Default-References)** – Information about the application steps when using inputted expressed data and default (ontology) references
  
* **3.2. [Input Expression Data and Uploaded References](https://github.com/AndorfLab/CytoPheno/wiki/3.2.-Input-Expression-Data-and-Uploaded-References)** – Information about the application steps when using inputted expressed data and uploaded (CSV file) references
  
**4. Input Marker Descriptors Tab** – Use when uploading or typing marker definitions (e.g. CD4+, CD8-) (Parts 2-3)  

* **4.1. [Input Marker Descriptors and Default References](https://github.com/AndorfLab/CytoPheno/wiki/4.1.-Input-Marker-Descriptors-and-Default-References)** – Information about the application steps when using inputted marker descriptions and default (ontology) references  

* **4.2. [Input Marker Descriptors and Uploaded References](https://github.com/AndorfLab/CytoPheno/wiki/4.2.-Input-Marker-Descriptors-and-Uploaded-References)** – Information about the application steps when using inputted marker descriptions and uploaded (CSV file) references  

**5. [Side-Panel Uploads and Parameters](https://github.com/AndorfLab/CytoPheno/wiki/5.-Side-Panel-Uploads-and-Parameters)** – Detailed info about all upload boxes and parameter widgets within the left-sided panel  

## Citation
[Tursi, A. R., Lages, C. S., Quayle, K., Koenig, Z. T., Loni, R., Eswar, S., Cobena-Reyes, J., Thornton, S., Tilburgs, T., & Andorf, S. Automated descriptive cell type naming in flow and mass cytometry with CytoPheno. Sci Rep 15, 26750 (2025). https://doi.org/10.1038/s41598-025-12153-w](https://doi.org/10.1038/s41598-025-12153-w)


