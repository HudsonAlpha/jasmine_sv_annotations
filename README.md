# Jasmine SV Annotator 
## Description
This is a pipeline that uses Jasmine (https://github.com/mkirsche/Jasmine) structural variant merging to (1) create a consensus SV set for an individual based on multiple SV callers (such as pbsv and sniffles) and then (2) create a merge between the individual consensus set and various population-scale structural variant databases (such as gnomad, hgsvc2, etc.) and (3) transfer the annotations of interest (such as population allele frequency) to the individual consensus set. 

## Usage
Currently this tool is intended for internal use at the Greg Cooper Lab at the HudsonAlpha Institute for Biotechnology. At this time it has not been written or tested for deployment outside of the HudsonAlpha compute cluster and is in active development. Thus, it may not run as an integrated pipeline in its current state. We are providing this repository primarily for publication transparency at this time (see Hiatt 2024 below).

## Hiatt 2024
Portions of this tool (steps 2 and 3 in Description) were used for descriptive structural variant analysis in Hiatt 2024 in Genome Research (Long-read genome sequencing and variant reanalysis increase diagnostic yield in neurodevelopmental disorders
https://doi.org/10.1101/gr.279227.124). Notes and scripts for that use case may be found in the hiatt_2024 subfolder. For now, please cite that manuscript if you use or modify this tool.

## Contact
James Lawlor, jlawlor@hudsonalpha.org
