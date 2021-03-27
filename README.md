# PGx_Guidelines
The univariate analysis consists of looking at drug sensitivity signatures for IC50 and AAC. 

# Download datasets
All of the utilized datasets for PGx Guidelines experiments are publicly available in the PSet format via ORCESTRA platform:
```
https://www.orcestra.ca/pset/stats
```
Specifically, CTRPv2, GDSCv2, CCLE, gCSI, and GDSCv1 are required for the univariate analysis.

# Run univariate analysis
Each Rscript includes code to load required libraries and datasets. 

Simply run the following for:
- all [solid and non-solid] tissues:
```
Rscript biomarker_analysis_alltissues.R "$@"
```

- after excluding non-solid tissues:
```
Rscript biomarker_analysis_solidonly.R "$@"
```

- after excluding non-solid tissues and log transformed IC50 values:
```
Rscript biomarker_analysis_log.R "$@"
```

- after excluding non-solid tissues and truncated
```
Rscript biomarker_analysis_truncated.R "$@"
```

- after excluding non-solid tissues, truncated, and log transformed IC50 values:
```
Rscript biomarker_analysis_truncated_log.R "$@"
```

