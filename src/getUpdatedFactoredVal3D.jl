function getUpdatedFactoredVal_3D(HO::Bool,kappaSquaredLoc::Float64,T::Array{EikTtype,3},Done::Array{Bool,3},loc1::Int64,
									loc2::Int64, loc3::Int64, h::Array{Float64,1}, hinv::Array{Float64,1}, n::Array{Int64,1},
									src::Array{Int64,1}, T0loc::Float64, G0loc::Array{Float64,1}, OpLoc::Array{Int8,1},A::Array{Float64,1},B::Array{Float64,1})
# op: -1 for backward. 1 for forward.
# A and B are just containers for answers so we won't allocate tuples all the time....
if HO
	calcDirectionHO_3D(1, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	calcDirectionHO_3D(2, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	calcDirectionHO_3D(3, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
else
	calcDirectionFO_3D(1, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	calcDirectionFO_3D(2, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
	calcDirectionFO_3D(3, T, Done, loc1, loc2, loc3, h, hinv, n, src, T0loc, G0loc, OpLoc,A,B);
end


#### Quadratic Solve
@inbounds A1 = A[1]; A2 = A[2]; A3 = A[3]; B1 = B[1]; B2 = B[2]; B3 = B[3];
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

ChooseBWD = getPaddedArr3D(Done,loc1m1,loc2m1,loc3m1);
ChooseFWD = getPaddedArr3D(Done,loc1p1,loc2p1,loc3p1);
if ChooseBWD && ChooseFWD
	@inbounds ChooseFWD = T[loc1p1,loc2p1,loc3p1]*T03D(n,h,src,loc1p1,loc2p1,loc3p1) < T[loc1m1,loc2m1,loc3m1]*T03D(n,h,src,loc1m1,loc2m1,loc3m1);
	ChooseBWD = !ChooseFWD;
end

A[dim] = 0.0;
B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD
	@inbounds if getPaddedArr3D(Done,loc1m2,loc2m2,loc3m2) #&& T[loc1m1,loc2m1,loc3m1]*T0[loc1m1,loc2m1,loc3m1] > T[loc1m2,loc2m2,loc3m2]*T0[loc1m2,loc2m2,loc3m2]
		@inbounds A[dim] = 1.5*t + G0loc[dim]; #going for second order!!!
        @inbounds B[dim] = (2*T[loc1m1,loc2m1,loc3m1] - 0.5*T[loc1m2,loc2m2,loc3m2])*t;
        @inbounds OpLoc[dim] = -2;
	else # first order.
		@inbounds A[dim] = t + G0loc[dim];
		@inbounds B[dim] = T[loc1m1,loc2m1,loc3m1]*t;
		@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
	end
elseif ChooseFWD
	@inbounds if getPaddedArr3D(Done,loc1p2,loc2p2,loc3p2) #&& T[loc1p1,loc2p1,loc3p1]*T0[loc1p1,loc2p1,loc3p1] > T[loc1p2,loc2p2,loc3p2]*T0[loc1p2,loc2p2,loc3p2]
		@inbounds A[dim] = 1.5*t - G0loc[dim];
        @inbounds B[dim] = (2*T[loc1p1,loc2p1,loc3p1] - 0.5*T[loc1p2,loc2p2,loc3p2])*t;
        @inbounds OpLoc[dim] = 2;
	else
		@inbounds A[dim] = t - G0loc[dim];
		@inbounds B[dim] = T[loc1p1,loc2p1,loc3p1]*t;
		@inbounds OpLoc[dim] = 1; 
	end
end
return;
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
A[dim] = 0.0;
B[dim] = 0.0;
OpLoc[dim] = 0;
t = T0loc*hinv[dim];
if ChooseBWD
	@inbounds A[dim] = t + G0loc[dim];
	@inbounds B[dim] = T[loc1m1,loc2m1,loc3m1]*t;
	@inbounds OpLoc[dim] = -1; # Choosing N means backward derivative.
elseif ChooseFWD
	@inbounds A[dim] = t - G0loc[dim];
	@inbounds B[dim] = T[loc1p1,loc2p1,loc3p1]*t;
	@inbounds OpLoc[dim] = 1; 
end

return;
end