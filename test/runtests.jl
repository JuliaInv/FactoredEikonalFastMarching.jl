using Test
using LinearAlgebra
using jInv.Mesh;
using FactoredEikonalFastMarching;
using Statistics
include("sensitivityTest.jl");

@testset "Sensitivity test" begin
	runSensitivityTest2D(); 
	runSensitivityTest3D();
	include("../examples/runExperiments.jl")
end
