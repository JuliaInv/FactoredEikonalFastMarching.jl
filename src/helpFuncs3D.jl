# 
# include("helpFuncs3D.jl"); n = [10,11,12]; loc1 = 3; loc2 = 4; loc3 = 5; (loc1a,loc2a,loc3a) = cs2loc3D(loc2cs3D(loc1,loc2,loc3,n),n)
# 

export loc2cs3D,cs2loc3D

function loc2cs3D(loc::Array{Int64,1},n::Array{Int64,1})
@inbounds cs = loc[1] + (loc[2]-1)*n[1] + (loc[3]-1)*n[1]*n[2];
return cs;
end

function loc2cs3D(loc1::Int64,loc2::Int64,loc3::Int64,n::Array{Int64,1})
@inbounds cs = loc1 + (loc2-1)*n[1] + (loc3-1)*n[1]*n[2];
return cs;
end

function cs2loc3D(loc::Array{Int64,1},cs_loc::Int64,n::Array{Int64,1})
@inbounds loc[1] = mod(cs_loc-1,n[1]) + 1;
@inbounds loc[2] = div(mod(cs_loc-1,n[1]*n[2]),n[1]) + 1;
@inbounds loc[3] = div(cs_loc-1,n[1]*n[2])+1;
return;
end

function cs2loc3D(cs_loc::Int64,n::Array{Int64,1})
@inbounds loc1 = mod(cs_loc-1,n[1])+1;
@inbounds loc2 = div(mod(cs_loc-1,n[1]*n[2]),n[1]) + 1;
@inbounds loc3 = div(cs_loc-1,n[1]*n[2])+1;
return (loc1,loc2,loc3);
end


function get3DNeighbors(loc::Array{Int64,1},n::Array{Int64,1},ans::Array{Int64,2},relevant::Array{Bool,1})
# This function fills the neighbours array ans with neighbours the "C" style, so it doesn't allocate.
@inbounds relevant[1]= loc[1]>1;
@inbounds relevant[2]= loc[1]<n[1];
@inbounds relevant[3]= loc[2]>1;
@inbounds relevant[4]= loc[2]<n[2];
@inbounds relevant[5]= loc[3]>1;
@inbounds relevant[6]= loc[3]<n[3];
######################################################################################
@inbounds ans[1,1] = loc[1]-1;
@inbounds ans[1,2] = loc[2];
@inbounds ans[1,3] = loc[3];
@inbounds ans[2,1] = loc[1]+1;
@inbounds ans[2,2] = loc[2];
@inbounds ans[2,3] = loc[3];
@inbounds ans[3,1] = loc[1];
@inbounds ans[3,2] = loc[2]-1;
@inbounds ans[3,3] = loc[3];
@inbounds ans[4,1] = loc[1];
@inbounds ans[4,2] = loc[2]+1;
@inbounds ans[4,3] = loc[3];
@inbounds ans[5,1] = loc[1];
@inbounds ans[5,2] = loc[2];
@inbounds ans[5,3] = loc[3]-1;
@inbounds ans[6,1] = loc[1];
@inbounds ans[6,2] = loc[2];
@inbounds ans[6,3] = loc[3]+1;
return;
end

function setPaddedArr3D(arr::Array{Bool,3},loc1::Int64,loc2::Int64,loc3::Int64,val::Bool)
# this function assumes that the boundaries are OK and the value is placed in a legal location.
@inbounds arr[loc1+2,loc2+2,loc3+2] = val
end

function getPaddedArr3D(arr::Array{Bool,3},loc1::Int64,loc2::Int64,loc3::Int64)
# this function assumes that the boundaries are OK and the value is placed in a legal location.
@inbounds return arr[loc1+2,loc2+2,loc3+2]
end
