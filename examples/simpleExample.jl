##
#
# This is a simple example for travel time calculation in 3D using the Fast Marching approach for the factored eikonal equation on a regular grid.
# See more details in:
# Eran Treister and Eldad Haber, A fast marching algorithm for the factored eikonal equation, Journal of Computational Physics, 324, 210-225, 2016.
#
# Please cite this paper if you use our package.
#
##
using jInv.Mesh
using FactoredEikonalFastMarching

plotting = false; 
if plotting
	using PyPlot;
	using jInvVis;
end


#### SETTING THE PARAMETERS OF THE PROBLEM

dim = 3;
n_nodes = 0;
if dim == 2
	n_nodes = [33;33];
	domain = [0.0,1.0,0.0,1.0]; # Domain is in KMs
	src = div(n_nodes,2);
	src[2] = 1;
else
	n_nodes = [33;33;33];
	domain = [0.0,1.0,0.0,1.0,0.0,1.0]; # Domain is in KMs
	src = div(n_nodes,2);
	src[3] = 1;
end
n_tup = tuple(n_nodes...);
mNodal = ones(Float64,n_tup); # m is 1/v^2 whre v is the velocity in km/sec
 
#### SETTING THE PARAMETERS OF THE FAST MARCHING SOLVER

M = getRegularMesh(domain,n_nodes-1);
mem = getEikonalTempMemory(n_nodes); # here we allocate temporary memory. Reuse this memory for multiple sources on the same mesh.
eikParam = getEikonalParam(M,mNodal,src,true); 
solveFastMarchingUpwindGrad(eikParam, mem);
mem = 0; # deallocate memory.
T = copy(eikParam.T1); # here we get the factor tau1.

## Below we extract the true travel time (multiply it wit tau0)
# This function works on the solution vector itself.
# Use getAnalytic2DeikonalSolutionAll and getAnalytic3DeikonalSolutionAll to get the analytical solution itself for a constant velocity model.
if M.dim==2 
	selfMultiplyWithAnalyticSolution2D(n_nodes,M.h,src,T);
else
	selfMultiplyWithAnalyticSolution3D(n_nodes,M.h,src,T)
end
# Here, T is the solution of the eikonal equation, the answer is in seconds.

#### PLOTTING THE RESULTS 

if plotting
	plotModel(T);
end

