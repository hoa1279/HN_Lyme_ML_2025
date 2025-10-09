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
1. Prepare fasta files including all aligned protein sequences and metadata file with classification of each sequence.
2. OHE_code.R: perform one-hot-encoding for protein sequences and save to "{gene}_OHE.csv".
3. Create_input_data.R: Prepare data for machine learning. Outputs are saved to "{gene}_fs_ML_input_data.RData".
4. ML_remote.R: Perform machine learning modelling.
5. Collect_model_performance.R: Collect all performance metrics, VIP scores, Top20 features.




------------------------------------------------------------------------

**Last Updated:** October 09, 2025\
**Version:** 1.0.0
