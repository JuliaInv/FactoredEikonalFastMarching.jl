function getUpdatedFactoredVal_2D(HO::Bool,kappa_squared::Array{Float64,2},T::Array{EikTtype,2},Done::Array{Bool,2},loc1::Int64,
									loc2::Int64,h::Array{Float64,1},hinv::Array{Float64,1},src::Array{Int64,1}, 
									T0loc::Float64, G0loc::Array{Float64,1}, OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
# op: -1 for backward. 1 for forward. -2/2 for 2nd order backward/forward. 
# A and B are just containers for answers so we won't allocate tuples all the time....
if HO
	#answers are in A,B,OpLoc
	(t1,t1a) = calcDirectionHO_2D(1, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
	(t2,t2a) = calcDirectionHO_2D(2, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
else
	#answers are in A,B,OpLoc
	t1 = calcDirectionFO_2D(1, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
	t2 = calcDirectionFO_2D(2, kappa_squared, T, Done, loc1, loc2, h, hinv, src, T0loc, G0loc, OpLoc, A, B);
end	

#### Quadratic Solve
A1 = A[1]; A2 = A[2]; B1 = B[1]; B2 = B[2];
a = A1*A1 + A2*A2; b = -A1*B1-A2*B2; c = B1*B1 + B2*B2; # b is divided by 2 here for easing the quadratic solve.
T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,kappa_squared[loc1,loc2],a,b,c,OpLoc);

if FactoredEikonalFastMarching.imposeMonotonicity
	if (t1 > 0.0 && T_loc*T0loc < t1+1e-16)
		if OpLoc[1]==-1 || OpLoc[1]==1
			A1 = T0loc/h[1];
			B1 = t1/h[1];	
		else
			A1 = (3.0*T0loc)/(2.0*h[1]);
			B1 = (4.0*t1 - t1a)/(2.0*h[1]);
		end
		a = A1*A1 + A2*A2; b = -A1*B1-A2*B2; c = B1*B1 + B2*B2; # b is divided by 2 here for easing the quadratic solve.
		T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,kappa_squared[loc1,loc2],a,b,c,OpLoc);
	elseif (t2 > 0.0 && T_loc*T0loc < t2+1e-16)
		if OpLoc[2]==-1 || OpLoc[2]==1
			A2 = T0loc/h[2];
			B2 = t2/h[2];
		else
			A2 = (3.0*T0loc)/(2.0*h[2]);
			B2 = (4.0*t2 - t2a)/(2.0*h[2]);
		end
		a = A1*A1 + A2*A2; b = -A1*B1-A2*B2; c = B1*B1 + B2*B2; # b is divided by 2 here for easing the quadratic solve.
		T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,kappa_squared[loc1,loc2],a,b,c,OpLoc);
	end
end 
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

tBWD = 0.0;
tFWD = 0.0;


ChooseBWD = getPaddedArr2D(Done,loc1m1,loc2m1);
if ChooseBWD
	tBWD = T[loc1m1,loc2m1]*T02D(h,src,loc1m1,loc2m1);
end
ChooseFWD = getPaddedArr2D(Done,loc1p1,loc2p1);
if ChooseFWD
	tFWD = T[loc1p1,loc2p1]*T02D(h,src,loc1p1,loc2p1)
end
if ChooseBWD && ChooseFWD
	ChooseFWD = tFWD < tBWD;
	ChooseBWD = !ChooseFWD;
end

A[dim] = B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD
	if getPaddedArr2D(Done,loc1m2,loc2m2) && tBWD > T[loc1m2,loc2m2]*T02D(h,src,loc1m2,loc2m2)
		@inbounds A[dim] = 1.5*t + G0loc[dim]; #going for second order!!!
        @inbounds B[dim] = (2.0*T[loc1m1,loc2m1] - 0.5*T[loc1m2,loc2m2])*t;
        @inbounds OpLoc[dim] = -2;
		return tBWD,T[loc1m2,loc2m2]*T02D(h,src,loc1m2,loc2m2);
	else # first order.
		@inbounds A[dim] = t + G0loc[dim];
		@inbounds B[dim] = T[loc1m1,loc2m1]*t;
		@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
		return tBWD,0.0;
	end
	
elseif ChooseFWD
	@inbounds if getPaddedArr2D(Done,loc1p2,loc2p2) && tFWD > T[loc1p2,loc2p2]*T02D(h,src,loc1p2,loc2p2)
		@inbounds A[dim] = 1.5*t - G0loc[dim];
        @inbounds B[dim] = (2.0*T[loc1p1,loc2p1] - 0.5*T[loc1p2,loc2p2])*t;
        @inbounds OpLoc[dim] = 2;
		return tFWD,T[loc1p2,loc2p2]*T02D(h,src,loc1p2,loc2p2);
	else
		@inbounds A[dim] = t - G0loc[dim];
		@inbounds B[dim] = T[loc1p1,loc2p1]*t;
		@inbounds OpLoc[dim] = 1;
		return tFWD,0.0;
	end
	
end
return 0.0,0.0;
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

tBWD = 0.0;
tFWD = 0.0;

ChooseBWD = getPaddedArr2D(Done,loc1m1,loc2m1);
ChooseFWD = getPaddedArr2D(Done,loc1p1,loc2p1);


if ChooseBWD && ChooseFWD
	tBWD = T[loc1m1,loc2m1]*T02D(h,src,loc1m1,loc2m1);
	tFWD = T[loc1p1,loc2p1]*T02D(h,src,loc1p1,loc2p1)
	@inbounds ChooseFWD = tFWD < tBWD;
	ChooseBWD = !ChooseFWD;
elseif FactoredEikonalFastMarching.imposeMonotonicity
	if ChooseBWD
		tBWD = T[loc1m1,loc2m1]*T02D(h,src,loc1m1,loc2m1);
	end
	if ChooseFWD
		tFWD = T[loc1p1,loc2p1]*T02D(h,src,loc1p1,loc2p1)
	end
end

@inbounds A[dim] = B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD		
	@inbounds A[dim] = t + G0loc[dim];
	@inbounds B[dim] = T[loc1m1,loc2m1]*t;
	@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
	return tBWD;
elseif ChooseFWD
	@inbounds A[dim] = t - G0loc[dim];
	@inbounds B[dim] = T[loc1p1,loc2p1]*t;
	@inbounds OpLoc[dim] = 1; 
	return tFWD;
end
return 0.0;
end