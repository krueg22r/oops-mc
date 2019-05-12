# oops-mc
Object-oriented polymer simulation (OOPS) package -- tools for performing Monte Carlo molecular simulation in multiple ensembles.

The OOPS package is a small, efficient c++ package designed for running flexible molecular simulations. Simulations may be performed in
the canonical (NVT) ensemble or the grand canonical (muVT) ensemble. The force field may be assembled from a combination of pair 
potentials (which depend on the relative position of two arbitrary particles), bond potentials (which depend on the relative positions
of particles in the same molecule), and external potentials (which depend on the position of each atom in the simulation box). 
The provided potentials may be used to perform coarse-grained bead-spring polymer simulations in periodic simulation boxes or 
under confinement. Several polymer-specific MC moves are included to improve sampling. The basic algorithm and the configurational
bias algorithm used in grand canonical simulations are described in Frenkel and Smit's Understanding Molecular Simulation. The Wikipedia page on 
Monte Carlo method in statistical physics gives a nice overview: 
https://en.wikipedia.org/wiki/Monte_Carlo_method_in_statistical_physics. 

## Prerequisites

The Eigen library (http://eigen.tuxfamily.org/index.php?title=Main_Page)
is used for linear algebra operations. 
Packmol (http://m3g.iqm.unicamp.br/packmol/home.shtml) 
is used by the input generator to create initial simulation 
configurations. 

## Installing

A sample Makefile is provided. The program has been tested using g++-4.9. 

## Running the program

The main input file uses keywords to supply simulation parameters and the names of files where molecular topology and initial 
coordinates are supplied. The order of keywords does not matter, except for sections specifying force field parameters. 
See the SampleInputs directory for several examples. The input file is provided as standard output. 

For example, to run the sample grand canonical bulk simulation on the command 
line, enter: 
./mc <in.gcBulk >out.gcBulk &

## Authors

Rachel Krueger

