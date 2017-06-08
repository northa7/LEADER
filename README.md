# LEADER

Latitudes and Elongations of Asteroid Distributions Estimated Rapidly (LEADER) is a software package designed to be used with **MATLAB**. The code has been tested on MATLAB version R2014a, so the R2014a version (or any newer release) should work as intended. For a detailed, technical documentation of the software package, it is strongly recommended to read Nortunen & Kaasalainen (2017).

The main branch of the software package has been written to work with the Wide-field Infrared Survey Explorer (WISE) database. In case there is interest in it, we may release the branch intended to work with the Panoramic Survey Telescope and Rapid Response System (Pan-STARRS) database.


## Files

In addition to this readme file and the license file, this software package includes the following files:
* ast_comparison_WISE.m - a routine for comparing asteroid populations
* damit_model.m - reads vertex and face information from a DAMIT asteroid model (for synthetic simulations)
* KS_comparison.m - a statistical comparison for two (p, beta) distributions
* lcg_read_synth_WISE.m - read JD + geometries from a WISE data file (for synthetic simulations)
* lcg_read_WISE.m - read JD + geometries + intensity from a WISE data file
* leader_brightness_synth_WISE.m - compute total brightness and eta (for synthetic simulations)
* leader_ellipsoid.m - compute semiaxes a, b and c for a DAMIT asteroid model (for synthetic simulations)
* leader_invert.m - construct analytical basis functions and compute the joint (p, beta) distribution
* leader_main_WISE.m - the main file for the primary routine of analyzing the characteristics of an asteroid population
* leader_phasecorr.m - a phase correction file for reducing the error caused by large changes in geometry
* leader_plots.m - plot the computed p and beta distributions
* leader_postprocess_WISE.m - perform a visual deconvolution to reduce the error of the computed (p, beta) distribution
* leader_synth_main_WISE.m - the main file for the synthetic simulations
* transform_mat.m - a coordinate transformation from the ecliptic frame to the asteroid's own frame
* weightbar.m - draw the density plot of a single variable (in bar form)

We have also included optional, separate zip archives of
* a small subpopulation of the WISE database (a sample to get started with the software package)
* asteroid models from DAMIT (a sample that can be used with synthetic simulations)
The archives contain separate installation instructions for their files.


## Usage

Each application is started by running their main function in MATLAB environment.

* Primary routine: run the code *leader_main_WISE.m* to compute the joint (p, beta) distribution for a (WISE) population. To use another subpopulation, create a variable called *lcg_files*, which contains the directory paths to the desired data files, and preload it to MATLAB's workspace before running *leader_main_WISE.m*.

* Synthetic simulator: run the code leader_synth_main_WISE.m to compute the joint (p, beta) distribution for a synthetic population. You need a file called *asteroideja.txt* in your MATLAB work directory, with the file containing the paths to object files containing the data (vertices and faces) for synthetic asteroid models. By default, the geometries are read from the WISE database. To use other another subpopulation for the database, create a variable called *lcg_files*, which contains the directory paths to the desired data files, and preload it to MATLAB's workspace before running *leader_synth_main_WISE.m*.

* Distribution comparator: run the code ast_comparison_WISE with two .mat-type files as input files. The .mat files must contain a variable called *lcg_files*, which contains the directory paths to the data files of the desired population. The function will load the contents of the .mat files into workspace and compute the (p, beta) distributions for both populations, and compare the obtained distributions. For example, to compare two populations with their data path files saved in files *WISE01.mat* and *WISE02.mat*, use the following syntax: *ast_comparison_WISE('WISE01', 'WISE02')*.

**You may freely edit and fine-tune the software package for personal use. Please note that if you experiment with your own database, it is strongly recommended that you run the synthetic simulator first, to test the accuracy and stability of the algorithm when using the database.**


## License

This software is licensed under [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/legalcode).

If you use LEADER in your research, please cite:
Nortunen, H.; Kaasalainen, M.; *LEADER: fast estimates of asteroid shape elongation and spin latitude distributions from scarce photometry*, A&A.


## Contact

Bug reports, feature suggestions, questions and comments are welcome. Please send them to Hari Nortunen (hari.nortunen@tut.fi).
