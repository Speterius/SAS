# Stochastic practical Matlab code

## Overview

All the code used in the analysis done for the AE4304P practical assignment simulating the longitudinal response of an aircraft to stochastic vertical gust input.

## Structure and running

Run *main.m*. A prompt will come up whether all figures should be opened. It might take a while for the figures to load up.

Each code section in *main.m* reflects one of the tasks in the assignment reader.
 
The *functions* folder contains relevant outside functions constructed. Most of them are based on the example code provided in the folder *ML_ae4304* (lecture note examples).

Relevant *functions*:

- *generate_state_space.m*
- *extend_state_space.m*
- *time_simulation.m*
- *analytic_psd.m*
- *experi_psd.m*

The *plot_psds.m* and *plot_psds_subplots.m* scripts are there to keep the code organized and readable (lots of lines of code) and are run from *main.m*.

## Authors

* **Peter Seres** - [P.Seres@student.tudelft.nl](P.Seres@student.tudelft.nl)

## Acknowledgments

* Supervisor: 	Dr. ir. Daan M. Pool
* Lecturer: 	Dr. ir. Max Mulder
* University:	Techinal University of Delft - Faculty of Aerospace Engineering


