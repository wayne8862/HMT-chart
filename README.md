# Holm’s multiple testing chart (HMT-chart)

The HMT-chart is a new Phase II control chart that is based on on a modified Holm's step-down multiple testing procedure (Holm (1979)) which achieves two important goals at the same time: (1) it simultaneously monitors correlated variables of different types, while keeping the probability of false alarm under desirable level, and (2) when the process is determined to be out of control, it further provides, without any additional efforts, diagnostics to pinpoint which parameters are out of control.


# Get Started
## Requirements
  - R version :  >= 3.6.3 
  - R package : mvtnorm, readr

## Simulation Study
 1. GenCorrelatedVariables.r : to generate the correlated random variables.
 2. Holm-and-Bon-type1-a-005.r : to simulate false alarm probabilities for the HMT-chart and BA-chart.
 3. HotellingT-type1-a-005.r : to simulate false alarm probabilities for the Hotelling T2-chart.
 4. simulation-Xpc-power-Case-all.r : to simulate the powers and the probability of correct diagnostics (p.c.d.) for the HMT-chart, BA-chart, and Hotelling T2-chart under the seven out-of-control scenarios.

## Data Example
 1. The data file "taiwan-monthly.csv" is in the **Data** folder.
 2. air-quality-monthly-resp.r : to run the data example application.

# Reference
Holm, S. (1979). A simple sequentially rejective multiple testing procedure, Scandinavian Journal of Statistics, 6, 65-70.
