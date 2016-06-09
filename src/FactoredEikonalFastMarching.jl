module FactoredEikonalFastMarching
export EikonalParam
export EikonalTempMemory
export getEikonalParam,getEikonalTempMemory
export solveFastMarchingUpwindGrad

using jInv.Mesh

EikIdxType  = Int32
EikTtype	= Float64 # here we choose the floating type for the calculation
export EikIdxType
const imposeMonotonicity = false;

"""
	type EikonalParam
	
	Stores results of a Factored Eikonal forward problem
	
	| nabla (T0(x) x T1(x)) | = m(x)
	
	where T0 is the analytical solution for m=1 and T1 is the factor 
	to be computed.
	
	Fields:
	
	Mesh 
	kappaSquared::Array{Float64}    - squared slowness model
	src         ::Array{Int64,1}    - index of source location
	HO          ::Bool              - flag for higher-order discretization
	T1          ::Array             - factor of Eikonal
	ordering    ::Array{EikIdxType} - order of points as were determined by 
	                                  Fast Marching (used in sensitivity computation)
	OP          ::Array{Int8}       - type of forward/backward operator in each 
	                                  voxel (used in sensitivity computation)

    Constructor: getEikonalParam
"""

type EikonalParam
    Mesh        :: RegularMesh # Mesh
	kappaSquared:: Array{Float64} # Squared slowness model
	src         :: Array{Int64,1}
	HO          :: Bool
	T1 			:: Array
	ordering	:: Array{EikIdxType,1}
	OP		    :: Array{Int8};
end


"""
type EikonalTempMemory
This struct is used for having a reusable memory allocation for consecutive runs of eikonal solves.
Done:: An array of booleans for knowing which variables are known.
J	:: An array of indices for the heap.
V   :: An array of values for the heap.
"""

type EikonalTempMemory
	Done		:: Array{Bool};
	J  			:: Array{Int64,1};
	V			:: Array{Float64,1};
end


function getEikonalParam(Mesh:: RegularMesh,kappaSquared::Array{Float64},src:: Array{Int64,1},HO:: Bool)
	if prod(Mesh.n+1) > (2^31-1)
		error("Eikonal is using Int32 for storing the sensitivities and N is too large. Change to EikIdxType in Eikoonal.jl to Int64");
	end
	pEik = EikonalParam(Mesh,kappaSquared,src,HO,[],[],[]);
	return pEik;
end

function getEikonalTempMemory(n::Array{Int64,1})
ntup = tuple(n...);
dim = length(n);
N = prod(n);
mem = EikonalTempMemory(zeros(Bool,tuple(n+4...)),zeros(Int64,N),zeros(Float64,N));
return mem;
end

import jInv.Utils.clear!
function clear!(param::EikonalParam)
param.kappaSquared = [];
param.T1 = [];
param.ordering = [];
param.OP = [];
return param;
end

function solveFastMarchingUpwindGrad(pEik::EikonalParam, mem::EikonalTempMemory)
if minimum(pEik.kappaSquared) < 1e-16
	error("Fast Marching::kappaSquared is negative somewhere");
end
 
dim = pEik.Mesh.dim;
if dim==3
	factoredFastMarching3D(pEik,mem);
elseif dim==2
	factoredFastMarching(pEik,mem);
end
return;
end

include("heap.jl");
include("helpFuncs2D.jl");
include("helpFuncs3D.jl");
include("factoredFastMarching2D.jl");
include("factoredFastMarching3D.jl");
include("getUpdatedFactoredVal2D.jl");
include("getUpdatedFactoredVal3D.jl");
include("getAnalyticEikonalSolution2D.jl");
include("getAnalyticEikonalSolution3D.jl");
include("EikSensitivities.jl");

end






