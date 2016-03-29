function getUpdatedFactoredVal_2D(HO::Bool,kappa_squared::Array{Float64,2},T::Array{EikTtype,2},Done::Array{Bool,2},loc1::Int64,
									loc2::Int64,h::Array{Float64,1},hinv::Array{Float64,1},src::Array{Int64,1}, 
									T0loc::Float64, G0loc::Array{Float64,1}, OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
# op: -1 for backward. 1 for forward. -2/2 for 2nd order backward/forward. 
# A and B are just containers for answers so we won't allocate tuples all the time....
if HO
	#answers are in A,B,OpLoc
	calcDirectionHO_2D(1, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
	calcDirectionHO_2D(2, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
else
	#answers are in A,B,OpLoc
	calcDirectionFO_2D(1, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
	calcDirectionFO_2D(2, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
end	

#### Quadratic Solve
A1 = A[1]; A2 = A[2]; B1 = B[1]; B2 = B[2];
a = A1*A1 + A2*A2; b = -A1*B1-A2*B2; c = B1*B1 + B2*B2; # b is divided by 2 here for easing the quadratic solve.
T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,kappa_squared[loc1,loc2],a,b,c,OpLoc);

return T_loc;
end

function solveQuadratic(a::Float64,b::Float64,c::Float64,fsq::Float64)
# Solution of a quadratic equation based on:
# t = (-b/2 + sqrt(b^2/4 - ac))/a. 
# b in the arguments is actually b/2 of the real equation. 

c -= fsq;
in_sqrt = b*b - a*c;
if (in_sqrt > 0.0)
    T_loc = (-b + sqrt(in_sqrt))/(a);
else
	# if in_sqrt < 0.0
		# warn("Eikonal: Negative value found in sqrt: ",in_sqrt," Solution probably has an error, because medium slowness is probably not smooth enough.");
	# end
    T_loc = 0.0;
end
return T_loc;
end

function solvePiecewiseQuadratic(A1::Float64,B1::Float64,A2::Float64,B2::Float64,fsq::Float64,
							a::Float64,b::Float64,c::Float64,OpLoc::Array{Int8,1})
T_loc = solveQuadratic(a,b,c,fsq);
g1 = (A1*T_loc - B1);
g2 = (A2*T_loc - B2);
if T_loc == 0.0 || g1 < 0.0 || g2 < 0.0
	val1 = (OpLoc[1]==0) ? Inf : B1/A1; 
	val2 = (OpLoc[2]==0) ? Inf : B2/A2;
	if val1 < val2
		f = sqrt(fsq);
		T_loc = (f + B1) / A1;
		OpLoc[2] = 0;
	else
		f = sqrt(fsq);
		T_loc = (f + B2) / A2;
		OpLoc[1] = 0;
	end
end


return T_loc
end

function calcDirectionHO_2D(dim::Int64,kappa_squared::Array{Float64,2},T::Array{EikTtype,2},Done::Array{Bool,2},
							loc1::Int64,loc2::Int64,h::Array{Float64,1},hinv::Array{Float64,1},src::Array{Int64,1},
							T0loc::Float64,G0loc::Array{Float64,1},OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
if dim==1
	loc1p1 = loc1 + 1;
	loc1m1 = loc1 - 1;
	loc1p2 = loc1 + 2;
	loc1m2 = loc1 - 2;
	loc2p1 = loc2m1 = loc2p2 = loc2m2 = loc2;
else # dim is assumed to be 2
	loc2p1 = loc2 + 1;
	loc2m1 = loc2 - 1;
	loc2p2 = loc2 + 2;
	loc2m2 = loc2 - 2;
	loc1p1 = loc1m1 = loc1p2 = loc1m2 = loc1;
end

ChooseBWD = getPaddedArr2D(Done,loc1m1,loc2m1);
ChooseFWD = getPaddedArr2D(Done,loc1p1,loc2p1);
if ChooseBWD && ChooseFWD
	@inbounds ChooseFWD = T[loc1p1,loc2p1]*T02D(h,src,loc1p1,loc2p1) < T[loc1m1,loc2m1]*T02D(h,src,loc1m1,loc2m1);
	ChooseBWD = !ChooseFWD;
end
A[dim] = B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD
	if getPaddedArr2D(Done,loc1m2,loc2m2) #&& abs( (kappa_squared[loc1m2,loc2m2] + kappa_squared[loc1,loc2]) - 2*kappa_squared[loc1m1,loc2m1]) < 0.01
		@inbounds A[dim] = 1.5*t + G0loc[dim]; #going for second order!!!
        @inbounds B[dim] = (2.0*T[loc1m1,loc2m1] - 0.5*T[loc1m2,loc2m2])*t;
        @inbounds OpLoc[dim] = -2;
	else # first order.
		@inbounds A[dim] = t + G0loc[dim];
		@inbounds B[dim] = T[loc1m1,loc2m1]*t;
		@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
	end
elseif ChooseFWD
	@inbounds if getPaddedArr2D(Done,loc1p2,loc2p2) #&& abs( (kappa_squared[loc1p2,loc2p2] + kappa_squared[loc1,loc2]) - 2*kappa_squared[loc1p1,loc2p1]) < 0.01
		@inbounds A[dim] = 1.5*t - G0loc[dim];
        @inbounds B[dim] = (2.0*T[loc1p1,loc2p1] - 0.5*T[loc1p2,loc2p2])*t;
        @inbounds OpLoc[dim] = 2;
	else
		@inbounds A[dim] = t - G0loc[dim];
		@inbounds B[dim] = T[loc1p1,loc2p1]*t;
		@inbounds OpLoc[dim] = 1; 
	end
end
return;
end

function calcDirectionFO_2D(dim::Int64,kappa_squared::Array{Float64,2},T::Array{EikTtype,2},Done::Array{Bool,2},loc1::Int64,loc2::Int64,
							h::Array{Float64,1},hinv::Array{Float64,1},src::Array{Int64,1},T0loc::Float64,G0loc::Array{Float64,1},
							OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})

if dim==1
	loc1p1 = loc1 + 1;
	loc1m1 = loc1 - 1;
	loc2p1 = loc2m1 = loc2;
else # dim is assumed to be 2
	loc2p1 = loc2 + 1;
	loc2m1 = loc2 - 1;
	loc1p1 = loc1m1 = loc1;
end

ChooseBWD = getPaddedArr2D(Done,loc1m1,loc2m1);
ChooseFWD = getPaddedArr2D(Done,loc1p1,loc2p1);
if ChooseBWD && ChooseFWD
	@inbounds ChooseFWD = T[loc1p1,loc2p1]*T02D(h,src,loc1p1,loc2p1) < T[loc1m1,loc2m1]*T02D(h,src,loc1m1,loc2m1);
	ChooseBWD = !ChooseFWD;
end
@inbounds A[dim] = B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD		
	@inbounds A[dim] = t + G0loc[dim];
	@inbounds B[dim] = T[loc1m1,loc2m1]*t;
	@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
elseif ChooseFWD
	@inbounds A[dim] = t - G0loc[dim];
	@inbounds B[dim] = T[loc1p1,loc2p1]*t;
	@inbounds OpLoc[dim] = 1; 
end
return;
end