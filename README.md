# Single-cell analysis reveals mitochondrial metabolism-linked AT2 stemness in bronchopulmonary dysplasia

Analysis code archive, version 0.1.0. This is a partial code release; end-to-end reproducibility has not been validated.

## Analysis layout

0-Phe: mitochondrial metabolism gene collection; 1-SC: single-cell QC, annotation and analysis; 2-AUCell: gene-set activity; 3-subset-AT2-CytoTRACE: AT2 subset; 4-scPagwas: GWAS integration; 5-CellChat: communication; 6-SC-GSVA: pathway analysis; 7-hdWGCNA: co-expression; 8-ML+SHAP: machine learning; 9-scTenifoldKnk: virtual knockout; 10-Ver: external validation.

Run scripts interactively from their containing directory so relative input/output paths resolve. Review scripts before execution: they include package installation, network calls and environment-specific settings. No single entry point or environment lockfile was provided. DEPENDENCIES.txt is extracted from source and is not a tested installation specification. The original session_info.txt is present in the source data folder under 7-hdWGCNA/output and should accompany the eventual data deposit.

## Data sources identified in supplied files

GEO: GSE217488, GSE225881, GSE220135, GSE188944 and GSE108754. FinnGen file: finngen_R12_P16_RESPI_DISTR_NEWBO.gz. These are file-derived identifiers; their role and inclusion in the final manuscript must be confirmed. Obtain upstream data from their original providers under the applicable terms. Data files are not included in this code candidate.

## Known missing materials

- 8-ML+SHAP/2-SHAP/source/my_fun.R referenced by the model training script.
- CytoTRACE analysis script: outputs exist, but no R script was supplied in 3-subset-AT2-CytoTRACE/2-CytoTRACE.
- Rps5 knockout main script: only plotKO.R was supplied for this target.
- Additional saved objects referenced in scripts may be missing; this repository has not been run end-to-end.

## Publication metadata

Authors (manuscript order): Rui Zhang, Jiaying Liu, Yu Gu, Li Fu, Bo Huang, and Zhenzhuang Zou.

Repository: https://github.com/baozizhuang/BPD-AT2-mitochondrial-stemness

Version: 0.1.0. Zenodo DOI: https://doi.org/10.5281/zenodo.24359512 (published code archive).

Original code is distributed under the MIT License; see LICENSE. Third-party code and materials retain their original terms, and this license does not override those terms. Data files are not included or licensed by this code archive. Credential-related lines were omitted from the publication copies; original source files are unchanged.
