# Comparing orthologies analysis

<p align="center"><img src="WorkflowCode.png" width="1000" /></p>


Manuscript: "Comparing orthology methods and their performance by recapitulating patterns of eukaryotic genome evolution."
Authors: E.S. Deutekom, B. Snel, T.J.P. van Dam

## Disclaimer
This serves to share code and data for reproducibility.
Due to the amount of data, only the files in green in the flowchart are provided. The flowchart shows a summary of the workflow, with the most important scripts and files. Some scripts that are used for supplementary figures and/or results are not shown. 
Any additional data can be shared upon request.

## Contained directories and files

Manuscript figures: Manuscript_Plots.ipynb

Supplementary figures: Supplementary Material Figures.ipynb

### eukarya directory
Contains part of the data and scripts used in the analyses minimally needed to reproduce the figures.

#### Software and tools in analyses
- blastp		2.7.1
- hhmbuild 		HMMER 3.1b2
- python		3.7.6
- samtools 		0.1.19-96b5f2294a

##### Python packages
- ete3			3.1.1
- biopython		1.76   
- numpy			1.18.1
- pandas		1.0.1
- scipy			1.4.1


### Figures directory
Contains the (supplementary) figures prduced for the manuscript.
Figures were produced using the following list of software and tools.

#### Software and tools
- python		3.7.6

##### Python packages
- matplotlib 		3.1.0
- pandas		0.24.2
- numpy			1.16.4
- scipy			1.3.0
- seaborn		0.10.0
- snakemake		5.5.4
