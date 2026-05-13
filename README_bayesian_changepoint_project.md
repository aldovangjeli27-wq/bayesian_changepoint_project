# Bayesian Changepoint Project

This repository contains an undergraduate Bayesian statistics project on DNA isochores and changepoint detection. The assignment studies whether each dataset is better described by **one** homogeneous isochore or by **two** isochores separated by a changepoint [file:51].

## Project idea

In the assignment, each dataset contains 100 consecutive observations. Every observation is the number of C or G bases inside a window of 5000 DNA bases, and the statistical question is whether the sequence comes from a single region with a stable CG rate or from two different regions with a structural change in that rate [file:51].

This repository approaches that question with Bayesian modeling. The code compares a one-isochore model against a two-isochore changepoint model and examines the posterior distribution of the changepoint location [file:48][file:51].

## Statistical setup

A Binomial model is used for the number of C/G bases in each 5000-base window, with Beta priors for the underlying probabilities. Under the one-isochore model, all windows share the same success probability; under the two-isochore model, the sequence is split at an unknown changepoint and each side has its own probability parameter [file:47][file:48][file:51].

The script `data_project.R` computes posterior quantities directly for the changepoint model, including posterior probabilities over candidate changepoint locations and model comparison quantities for one versus two isochores [file:48].

The script `R-project-gibbs.R` implements a Gibbs sampler for the two-isochore model. It samples the changepoint and the two Binomial probabilities, while using log-scale normalization and probability clipping for numerical stability [file:47].

## What makes this project interesting

A practical difficulty in this assignment was that the dataset could lead the Gibbs sampler to poor initial regions and unstable mixing behavior. To address this, the implementation includes an initialization heuristic (`pick_start_k`) that searches over candidate starting points and favors higher-posterior regions before running the main sampler [file:47].

This was developed as a debugging solution to a real sampling problem and later recognized as being related to standard collapsed or partially collapsed ideas used to improve MCMC behavior. The value of the project is not only the final model, but also the attention paid to convergence, initialization, and numerical robustness [file:47].

## Repository contents

- `data_project.R`: direct posterior analysis for the changepoint problem, including posterior probabilities for candidate changepoints [file:48]
- `R-project-gibbs.R`: Gibbs sampler implementation for the two-isochore model [file:47]
- `README.md`: overview of the project and modeling choices

## Main goals

- Model the number of C/G bases per window using a Bayesian Binomial-Beta framework [file:51]
- Compare a one-isochore model with a two-isochore changepoint model [file:51]
- Estimate the posterior distribution of the changepoint location [file:48]
- Explore practical Gibbs sampling issues such as initialization and mixing [file:47]

## Notes

This repository reflects coursework from Bayesian Statistics and Applications. The project was designed to answer a concrete scientific question about DNA structure while also serving as a hands-on exercise in Bayesian model comparison and MCMC implementation [file:51].

Some parts of the code are intentionally kept close to the original coursework version. Future improvements could include cleaner function structure, additional diagnostics, trace plots, and a more polished presentation of results [file:47][file:48].
