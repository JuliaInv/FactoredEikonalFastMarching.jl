
export loc2cs,cs2loc

function loc2cs(loc::Array{Int64,1},n::Array{Int64,1})
@inbounds cs = loc[1] + (loc[2]-1)*n[1];
return cs;
end

function loc2cs(loc1::Int64,loc2::Int64,n::Array{Int64,1})
@inbounds cs = loc1 + (loc2-1)*n[1];
return cs;
end

function cs2loc(loc::Array{Int64,1},cs_loc::Int64,n::Array{Int64,1})
@inbounds loc[1] = mod(cs_loc-1,n[1]) + 1;
@inbounds loc[2] = div(cs_loc-1,n[1])+1;
return;
end

function cs2loc(cs_loc::Int64,n::Array{Int64,1})
@inbounds loc1 = mod(cs_loc-1,n[1]) + 1;
@inbounds loc2 = div(cs_loc-1,n[1])+1; 
return loc1,loc2;
end


function get2DNeighbors(loc::Array{Int64,1},n::Array{Int64,1},ans::Array{Int64,2},relevant::Array{Bool,1})
# This function fills the neighbours array ans with neighbours the "C" style, so it doesn't allocate.
@inbounds relevant[1]= loc[1]>1;
@inbounds relevant[2]= loc[1]<n[1];
@inbounds relevant[3]= loc[2]>1;
@inbounds relevant[4]= loc[2]<n[2];

@inbounds ans[1,1] = loc[1]-1;
@inbounds ans[1,2] = loc[2];
@inbounds ans[2,1] = loc[1]+1;
@inbounds ans[2,2] = loc[2];
@inbounds ans[3,1] = loc[1];
@inbounds ans[3,2] = loc[2]-1;
@inbounds ans[4,1] = loc[1];
@inbounds ans[4,2] = loc[2]+1;
return;
end



function setPaddedArr2D(arr::Array{Bool,2},loc1::Int64,loc2::Int64,val::Bool)
# this function assumes that the boundaries are OK and the value is placed in a legal location.
@inbounds arr[loc1+2,loc2+2] = val
end


function getPaddedArr2D(arr::Array{Bool,2},loc1::Int64,loc2::Int64)
# this function assumes that the boundaries are OK and the value is placed in a legal location.
@inbounds return arr[loc1+2,loc2+2]
end

function ddxLongDiff(h::Float64,n::Int64)
	dm1 = (-1/(2*h)*ones(n-1));
	dm1[n-1] = -1/h;
	d   = zeros(n);
	d[1] = -1/h;
	d[n] = 1/h;
	dp1 = (1/(2*h))*ones(n-1);
	dp1[1] = 1/h;
	D = spdiagm((dm1,d,dp1),(-1,0,1),n,n);
	return D;
end

