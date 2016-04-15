function getUpdatedFactoredVal_3D(HO::Bool,kappaSquaredLoc::Float64,T::Array{EikTtype,3},Done::Array{Bool,3},loc1::Int64,
									loc2::Int64, loc3::Int64, h::Array{Float64,1}, hinv::Array{Float64,1}, n::Array{Int64,1},
									src::Array{Int64,1}, T0loc::Float64, G0loc::Array{Float64,1}, OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
# op: -1 for backward. 1 for forward.
# A and B are just containers for answers so we won't allocate tuples all the time....
if HO
	(t1,t1a) = calcDirectionHO_3D(1, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	(t2,t2a) = calcDirectionHO_3D(2, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	(t3,t3a) = calcDirectionHO_3D(3, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
else
	t1 = calcDirectionFO_3D(1, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	t2 = calcDirectionFO_3D(2, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	t3 = calcDirectionFO_3D(3, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
end

#### Piecewise Quadratic Solve
@inbounds A1 = A[1]; A2 = A[2]; A3 = A[3]; B1 = B[1]; B2 = B[2]; B3 = B[3];
T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,A3,B3,kappaSquaredLoc,OpLoc);


if FactoredEikonalFastMarching.imposeMonotonicity
	if (t1 > 0.0 && T_loc*T0loc < t1+1e-16)
		if OpLoc[1]==-1 || OpLoc[1]==1
			A1 = T0loc/h[1];
			B1 = t1/h[1];	
		else
			A1 = (3.0*T0loc)/(2.0*h[1]);
			B1 = (4.0*t1 - t1a)/(2.0*h[1]);
		end
		# A1 = 0.0;
		# B1 = 0.0;
		# OpLoc[1] = 0;
	end
	if (t2 > 0.0 && T_loc*T0loc < t2+1e-16)
		if OpLoc[2]==-1 || OpLoc[2]==1
			A2 = T0loc/h[2];
			B2 = t2/h[2];
		else
			A2 = (3.0*T0loc)/(2.0*h[2]);
			B2 = (4.0*t2 - t2a)/(2.0*h[2]);
		end
		# A2 = 0.0;
		# B2 = 0.0;
		# OpLoc[2] = 0;
	end
	if (t3 > 0.0 && T_loc*T0loc < t3+1e-16)
		if OpLoc[3]==-1 || OpLoc[3]==1
			A3 = T0loc/h[3];
			B3 = t3/h[3];
		else
			A3 = (3.0*T0loc)/(2.0*h[3]);
			B3 = (4.0*t3 - t3a)/(2.0*h[3]);
		end
		#A3 = 0.0;
		#B3 = 0.0;
		#OpLoc[3] = 0;
	end
	T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,A3,B3,kappaSquaredLoc,OpLoc);
end
return T_loc;
end


function solvePiecewiseQuadratic(A1::Float64,B1::Float64,A2::Float64,B2::Float64,A3::Float64,B3::Float64,kappaSquaredLoc::Float64,OpLoc::Array{Int8,1})
a1 = A1*A1; a2 = A2*A2; a3 = A3*A3; b1 = -A1*B1; b2 = -A2*B2; b3 = -A3*B3; c1 = B1*B1; c2 = B2*B2; c3 = B3*B3;
T_loc = solveQuadratic(a1+a2+a3,b1+b2+b3,c1+c2+c3,kappaSquaredLoc);


g1 = (A1*T_loc - B1);
g2 = (A2*T_loc - B2);
g3 = (A3*T_loc - B3);


if g1<0.0 || g2<0.0 || g3<0.0
	val1 = (OpLoc[1]==0) ? Inf : B1/A1; 
	val2 = (OpLoc[2]==0) ? Inf : B2/A2;
	val3 = (OpLoc[3]==0) ? Inf : B3/A3;
	# Below, solvePiecewiseQuadratic uses only the first two cells in OpLoc.
	if val1 < val2 # val1 is in for two.
		if val2 < val3 # solve for val1 and val2.
			T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,kappaSquaredLoc,a1+a2,b1+b2,c1+c2,OpLoc);
			OpLoc[3] = 0; #G_loc3 = 0.0;
		else # solve for val1 and val3.
			OpLoc[2] = OpLoc[3];
			T_loc = solvePiecewiseQuadratic(A1,B1,A3,B3,kappaSquaredLoc,a1+a3,b1+b3,c1+c3,OpLoc);
			OpLoc[3] = OpLoc[2]; # this is a hack to handle dimension 1 and 3 by a 2D code.
			OpLoc[2] = 0; #G_loc2 = 0.0;
		end
	else # val2 is in for two
		if val1 < val3 # solve for val2 and val1
 			T_loc = solvePiecewiseQuadratic(A1,B1,A2,B2,kappaSquaredLoc,a1+a2,b1+b2,c1+c2,OpLoc);
			OpLoc[3] = 0; #G_loc3 = 0.0;
		else # solve for val2 and val3.
		    OpLoc[1] = OpLoc[3];
			T_loc = solvePiecewiseQuadratic(A3,B3,A2,B2,kappaSquaredLoc,a3+a2,b3+b2,c3+c2,OpLoc);
			OpLoc[3] = OpLoc[1];
			OpLoc[1] = 0; #G_loc1 = 0.0;
		end
	end
end
return T_loc;
end


function calcDirectionHO_3D(dim::Int64,T::Array{EikTtype,3},Done::Array{Bool,3},
	loc1::Int64,loc2::Int64,loc3::Int64,h::Array{Float64,1},hinv::Array{Float64,1},n::Array{Int64,1},src::Array{Int64,1},
	T0loc::Float64,G0loc::Array{Float64,1},OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
# OpLoc is used as an output to avoid tuple allocation.
if dim == 1
	loc1p1 = loc1 + 1;
	loc1m1 = loc1 - 1;
	loc1p2 = loc1 + 2;
	loc1m2 = loc1 - 2;
	loc2p1 = loc2m1 = loc2p2 = loc2m2 = loc2;
	loc3p1 = loc3m1 = loc3p2 = loc3m2 = loc3;
elseif dim == 2 
	loc2p1 = loc2 + 1;
	loc2m1 = loc2 - 1;
	loc2p2 = loc2 + 2;
	loc2m2 = loc2 - 2;
	loc1p1 = loc1m1 = loc1p2 = loc1m2 = loc1;
	loc3p1 = loc3m1 = loc3p2 = loc3m2 = loc3;
else
	loc3p1 = loc3 + 1;
	loc3m1 = loc3 - 1;
	loc3p2 = loc3 + 2;
	loc3m2 = loc3 - 2;
	loc1p1 = loc1m1 = loc1p2 = loc1m2 = loc1;
	loc2p1 = loc2m1 = loc2p2 = loc2m2 = loc2;
end

tBWD = 0.0;
tFWD = 0.0;


ChooseBWD = getPaddedArr3D(Done,loc1m1,loc2m1,loc3m1);
if ChooseBWD
	tBWD = T[loc1m1,loc2m1,loc3m1]*T03D(n,h,src,loc1m1,loc2m1,loc3m1);
end
ChooseFWD = getPaddedArr3D(Done,loc1p1,loc2p1,loc3p1);
if ChooseFWD
	tFWD = T[loc1p1,loc2p1,loc3p1]*T03D(n,h,src,loc1p1,loc2p1,loc3p1);
end
if ChooseBWD && ChooseFWD
	ChooseFWD = tFWD < tBWD;
	ChooseBWD = !ChooseFWD;
end

A[dim] = 0.0;
B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD
	@inbounds if getPaddedArr3D(Done,loc1m2,loc2m2,loc3m2) && tBWD > T[loc1m2,loc2m2,loc3m2]*T03D(n,h,src,loc1m2,loc2m2,loc3m2)
		@inbounds A[dim] = 1.5*t + G0loc[dim]; #going for second order!!!
        @inbounds B[dim] = (2*T[loc1m1,loc2m1,loc3m1] - 0.5*T[loc1m2,loc2m2,loc3m2])*t;
        @inbounds OpLoc[dim] = -2;
		return (tBWD,T[loc1m2,loc2m2,loc3m2]*T03D(n,h,src,loc1m2,loc2m2,loc3m2))
	else # first order.
		@inbounds A[dim] = t + G0loc[dim];
		@inbounds B[dim] = T[loc1m1,loc2m1,loc3m1]*t;
		@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
		return (tBWD,0.0)
	end
elseif ChooseFWD
	@inbounds if getPaddedArr3D(Done,loc1p2,loc2p2,loc3p2) && tFWD > T[loc1p2,loc2p2,loc3p2]*T03D(n,h,src,loc1p2,loc2p2,loc3p2)
		@inbounds A[dim] = 1.5*t - G0loc[dim];
        @inbounds B[dim] = (2*T[loc1p1,loc2p1,loc3p1] - 0.5*T[loc1p2,loc2p2,loc3p2])*t;
        @inbounds OpLoc[dim] = 2;
		return (tFWD,T[loc1p2,loc2p2,loc3p2]*T03D(n,h,src,loc1p2,loc2p2,loc3p2));
	else
		@inbounds A[dim] = t - G0loc[dim];
		@inbounds B[dim] = T[loc1p1,loc2p1,loc3p1]*t;
		@inbounds OpLoc[dim] = 1; 
		return (tFWD,0.0);
	end
end
return (0.0,0.0);
end



function calcDirectionFO_3D(dim::Int64,T::Array{EikTtype,3},Done::Array{Bool,3},loc1::Int64,loc2::Int64,loc3::Int64,
h::Array{Float64,1},hinv::Array{Float64,1},n::Array{Int64,1},src::Array{Int64,1},T0loc::Float64,G0loc::Array{Float64,1},OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
loc1p1 = loc1m1 = loc1;
loc2p1 = loc2m1 = loc2;
loc3p1 = loc3m1 = loc3;
if dim == 1
	loc1p1 += 1;
	loc1m1 -= 1;
elseif dim == 2 # dim is assumed to be 2
	loc2p1 += 1;
	loc2m1 -= 1;
else
	loc3p1 += 1;
	loc3m1 -= 1;
end

ChooseBWD = getPaddedArr3D(Done,loc1m1,loc2m1,loc3m1);
ChooseFWD = getPaddedArr3D(Done,loc1p1,loc2p1,loc3p1);
if ChooseBWD && ChooseFWD
	@inbounds ChooseFWD = T03D(n,h,src,loc1p1,loc2p1,loc3p1)*T[loc1p1,loc2p1,loc3p1] < T[loc1m1,loc2m1,loc3m1]*T03D(n,h,src,loc1m1,loc2m1,loc3m1);
	ChooseBWD = !ChooseFWD;
end


tBWD = 0.0;
tFWD = 0.0;

ChooseBWD = getPaddedArr3D(Done,loc1m1,loc2m1,loc3m1);
ChooseFWD = getPaddedArr3D(Done,loc1p1,loc2p1,loc3p1);


if ChooseBWD && ChooseFWD
	tBWD = T[loc1m1,loc2m1,loc3m1]*T03D(n,h,src,loc1m1,loc2m1,loc3m1);
	tFWD = T03D(n,h,src,loc1p1,loc2p1,loc3p1)*T[loc1p1,loc2p1,loc3p1];
	@inbounds ChooseFWD = tFWD < tBWD;
	ChooseBWD = !ChooseFWD;
elseif FactoredEikonalFastMarching.imposeMonotonicity
	if ChooseBWD
		tBWD = T[loc1m1,loc2m1,loc3m1]*T03D(n,h,src,loc1m1,loc2m1,loc3m1);
	end
	if ChooseFWD
		tFWD =  T03D(n,h,src,loc1p1,loc2p1,loc3p1)*T[loc1p1,loc2p1,loc3p1];
	end
end

A[dim] = 0.0;
B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD
	@inbounds A[dim] = t + G0loc[dim];
	@inbounds B[dim] = T[loc1m1,loc2m1,loc3m1]*t;
	@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
	return tBWD;
elseif ChooseFWD
	@inbounds A[dim] = t - G0loc[dim];
	@inbounds B[dim] = T[loc1p1,loc2p1,loc3p1]*t;
	@inbounds OpLoc[dim] = 1; 
	return tFWD;
end
return 0.0;
end
