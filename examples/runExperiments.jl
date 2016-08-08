
using jInv.Mesh;
using FactoredEikonalFastMarching;
using PyPlot
close("all")

include("runAccuracyExperiment.jl");
include("getAnalyticalMediums.jl");
include("getWorkUnit.jl");

function runPaperExperiments()
numOfRefinements = 4; 

println("******************** 2D experiments ***************************");

I = [4,8];
n = zeros(Int64,2);

h0 = [0.05,0.05];
WU = zeros(numOfRefinements)
for k = 1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	WU[k] = GetWorkunitForLoop(n,h);
end

println("******************** Gradient Slowness ***************************");
println("h\t\tn\t\terror(1st)[max,l2]\ttime(1st)\terror(2nd)[max,l2]\ttime(2nd)")



h0 = [0.05,0.05];
for k=1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	(kappaSquared,src,T_exact) = getAnalyticalConstGrad2D(n,h);
	runExperimentAndWriteResults(kappaSquared,h,src,n,T_exact,WU[k]);
end

println("******************** Gradient Velocity ***************************");

println("h\t\tn\t\terror(1st)[max,l2]\ttime(1st)\terror(2nd)[max,l2]\ttime(2nd)")
h0 = [0.05,0.05];
for k=1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	(kappaSquared,src,T_exact) = getAnalyticalConstGradInv2D(n,h);
	runExperimentAndWriteResults(kappaSquared,h,src,n,T_exact,WU[k]);
end



println("******************** Gaussian factor ***************************");

println("h\t\tn\t\terror(1st)[max,l2]\ttime(1st)\terror(2nd)[max,l2]\ttime(2nd)")
h0 = [0.05,0.05];
for k=1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	(kappaSquared,src,T_exact) = getSmoothGaussianMedium(n,h);
	runExperimentAndWriteResults(kappaSquared,h,src,n,T_exact,WU[k]);
end

println("*******************************************************************************")
println("*******************************************************************************")
println("************************* 3-D EXPERIMENTS *************************************")
println("*******************************************************************************")
println("*******************************************************************************")

##########################################################################################
##########################################################################################
# 3D experiments:
##########################################################################################

I = [1.6,1.6,0.8];
n = zeros(Int64,3);

numOfRefinements = 4;

h0 = [0.1,0.1,0.1];
WU = zeros(numOfRefinements)
for k = 1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	n[3] = round(Int64,(I[3]/h[3])+1);
	WU[k] = GetWorkunitForLoop(n,h);
end

println("******************** Gradient Slowness ***************************");
println("h\t\tn\t\terror(1st)[max,l2]\ttime(1st)\terror(2nd)[max,l2]\ttime(2nd)")

for k=1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	n[3] = round(Int64,(I[3]/h[3])+1);
	(kappaSquared,src,T_exact) = getAnalyticalConstGrad3D(n,h);
	runExperimentAndWriteResults(kappaSquared,h,src,n,T_exact,WU[k]);
end

println("******************** Gradient Velocity ***************************");
println("h\t\tn\t\terror(1st)[max,l2]\ttime(1st)\terror(2nd)[max,l2]\ttime(2nd)")
for k=1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	n[3] = round(Int64,(I[3]/h[3])+1);
	(kappaSquared,src,T_exact) = getAnalyticalConstGradInv3D(n,h);
	runExperimentAndWriteResults(kappaSquared,h,src,n,T_exact,WU[k]);
end

println("******************** Gaussian factor ***************************");
println("h\t\tn\t\terror(1st)[max,l2]\ttime(1st)\terror(2nd)[max,l2]\ttime(2nd)")
for k=1:numOfRefinements
	h = (0.5^k)*h0;
	n[1] = round(Int64,(I[1]/h[1])+1);
	n[2] = round(Int64,(I[2]/h[2])+1);
	n[3] = round(Int64,(I[3]/h[3])+1);
	(kappaSquared,src,T_exact) = getSmoothGaussianMedium(n,h);
	runExperimentAndWriteResults(kappaSquared,h,src,n,T_exact,WU[k]);
end
return;
end
############################################################################################

runPaperExperiments();
