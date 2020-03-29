# maths-project

This repository contains the source code used to retrieve and preprocess global satellite wildfire data (predominantly in Python), as well as to perform statistical analysis on it (in R), for my Mathematical Project in partial fulfillment of my degree.

The `utils/` folder comprises quick and dirty Python scripts that I wrote to retrieve huge amounts of data (~200 GB) online.

The `preprocess.R` script processes various `.nc` files retrieved and turns it into a `.csv` file that is read by `model_fitting.R` to fit our Poisson hurdle model with latent spatial effects.

`model_checking.R` produces several model diagnostics and plots to examine the validity our model, and `simulation.R` performs a simulation study (not included in the report) using the parameters found in the model.

![Plot of daily average wildfire radiative power for Indonesia in 2015](/output/frpfire.ina.2015.gif)
