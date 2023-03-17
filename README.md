# GenerativeEmotionCodeArchive
Code archive to reproduce all analysis, simulations and figures in the paper "Towards a Generative Model for Emotion Dynamics".

# Overview
To run the code, install the `GenerativeModel` R package using:

```{r}
devtools::install_github('https://github.com/ryanoisin/GenerativeEmotion')
```

This repository includes the following files:

- `analysis_figures.R`: Creates all figures and stored them in `/figures`.
- `sensitivity_analysis.R`: Includes code or the sensitivity analysis. reported in the appendix of the paper.
- `aux_function.R`: Includes auxiliary function.
- `/files`: Includes data from Rowland & Wenzel (2020) used in `analysis_figures.R`.
