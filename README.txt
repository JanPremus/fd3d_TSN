FD3D_TSN code
 Authors: Jan Premus and Frantisek Gallovic (8/2019)
 Charles University in Prague, Faculty of Mathematics and Physics

 This code is published under the GNU General Public License. To any
 licensee is given permission to modify the work, as well as to copy
 and redistribute the work or any derivative version. Still we would
 like to kindly ask you to acknowledge the authors and don't remove
 their names from the code. This code is distributed in the hope
 that it will be useful, but WITHOUT ANY WARRANTY.
 ------------------------------------------------------

Dynamic simulation of spontaneous rupture propagation on planar fault using finite difference approach. 
The aim of the development is to provide fast and easy to use tool for earthquake source simulations.
Optional GPU accelearation is a part of the code.
To compile the code use INTEL FORTRAN compiler (https://software.intel.com/en-us/fortran-compilers) for the CPU, or
PGI FORTRAN compiler (https://www.pgroup.com/products/community.htm) for the GPU accelerated version.

Several examples of the use of the code are shown in the folder examples.

Amatrice 2016 earthquake - best fitting forward model from dynamic inversion (Gallovic et al., 2019).
	- Copy all source files into the folder Amatrice2016
	- Run compile.sh for compiling the code
	- Run fd3d_pt_CPU_TSN for the simulation

TPV5 benchmark - One of the community SCEC/USGS Spontaneous Rupture Code Verification Project (http://scecdata.usc.edu/cvws/)
benchmarks (Harris et al., 2019) - slip weakening friction law with heterogenous traction
	- Copy all source files into the folder TPV5
	- Run compile.sh for compiling the code (compiler flags -DSCEC and -DTPV5 are set in the file to properly set up and run the benchmark)
	- Run fd3d_pt_CPU_TSN for the simulation

TPV104 benchmark - One of the community SCEC/USGS Spontaneous Rupture Code Verification Project (http://scecdata.usc.edu/cvws/)
benchmarks (Harris et al., 2019) - fast velocity weakening friction law with heterogenous traction
	- Copy all source files into the folder TPV104 (compiler flags -DFWV -DSCEC and -DTPV104 are set in the file to properly set up the friction law and run the benchmark)
	- Run compile.sh for compiling the code
	- Run fd3d_pt_CPU_TSN for the simulation
TPV9 benchmark - One of the community SCEC/USGS Spontaneous Rupture Code Verification Project (http://scecdata.usc.edu/cvws/)
benchmarks (Harris et al., 2019) - fast velocity weakening friction law with heterogenous traction
	- Copy all source files into the folder TPV104 (compiler flags -DSCEC and -DTPV9 are set in the file to properly set up the friction law and run the benchmark)
	- Run compile.sh for compiling the code
	- Run fd3d_pt_CPU_TSN for the simulation

Commentary in compile.sh files contains a line to compile the GPU accelerated version. 
Switching between slip weakening (default) and fast velocity weakening friction is done using compiler flag -DFVW.
There are several outputs in the folder result:
	- sliprateX.res, sliprateZ.res - time histories of horizontal and vertical slip rate at all on fault nodes.
	- shearstressX.res, shearstressZ.res -  time histories of horizontal and vertical traction at all on fault nodes.
	- when seismograms are desired, positions of the seismic stations need to be provided in file inputfd3d.dat
		- set number of seismic stations (default 0) 
		- set positions of seismic stations on next lines (one line per station, integers): horizontal, normal and vertical position in FD grid, fault is at the normal position nyt
		- read seismograms in files stan$X$.res
	- Two MATLAB/OCTAVE files are provided to see the results
		- PrintSnapshot.m - prints graph of spatial distribution of slip rate and traction at a given time T
		- PrintTimeSeries.m - prints time series of slip rate and traction at a given node
Further documentation is provided in the doc folder. File FD3D_TSN Theoretical background describes the implementation of the finite difference method and boundary conditions, 
while FD3D_TSN Using the code describes input/output files and compilation options. 
