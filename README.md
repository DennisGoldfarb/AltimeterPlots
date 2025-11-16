# AltimeterPlots

This repository hosts a minimal Shiny application that queries the Koina Altimeter spline model to retrieve peptide fragmentation predictions. The current UI focuses on submitting a single peptide sequence and displaying the resulting prediction table. Future work can expand on these outputs to add interactive visualizations.

## Requirements

- R (>= 4.1 recommended)
- The `shiny` and `koinar` packages. Install via:

```r
install.packages("shiny")
remotes::install_github("wilhelm-lab/koinar")
```

## Running the app

From the repository root, run:

```r
shiny::runApp()
```

This will launch the application in your default browser. Provide a peptide sequence and a precursor charge (between 1 and 7) and click **Predict** to submit a request to Koina. The predictions will be displayed as a table until additional visualization layers are added.
