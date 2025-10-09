# Machine Learning-Driven Identification of Virulence Determinants in Borrelia burgdorferi Associated with Human Dissemination

## Overview

This repository contains the R code and analysis scripts for the
manuscript "Machine Learning-Driven Identification of Virulence
Determinants in Borrelia burgdorferi Associated with Human
Dissemination" by Hoa T. Nguyen and Catherine A. Brissette (2025).

**Paper Status:** \
**Journal:**\
**DOI:**\
**Preprint:** [10.1101/2025.07.09.663762](https://doi.org/10.1101/2025.07.09.663762)

## Abstract

This study applies machine learning (ML) approaches to predict Lyme
disease dissemination phenotypes from Borrelia burgdorferi protein
sequences. Using whole genome data from 299 clinical isolates, we
identified specific amino acid residues in surface-exposed virulence
factors that distinguish between localized and disseminated infections.
Our computational framework achieved robust predictive performance
(\>0.7 for all metrics) for three key proteins (DbpA, OspC, and RevA)
and revealed that ML-identified residues are significantly enriched in
predicted B-cell epitope regions, suggesting their role in immune
evasion and bacterial persistence

## Contact

For questions about this code or research:

-   **Corresponding Author:** Hoa T. Nguyen
    ([hoa.t.nguyen\@und.edu](mailto:hoa.t.nguyen@und.edu))



## Execution order:
1. Prepare FASTA files (aligned protein sequences) and metadata file (sequence classifications). Files can be found in `data/{gene}/` folder.
2. OHE_code.R - One-hot encode sequences → {gene}_OHE.csv
3. Create_input_data.R - Prepare ML input data → {gene}_fs_ML_input_data.RData
4. ML_remote.R - Train machine learning models
5. Collect_model_performance.R - Compile performance metrics, VIP scores, and top features




------------------------------------------------------------------------

**Last Updated:** October 09, 2025\
**Version:** 1.0.0
