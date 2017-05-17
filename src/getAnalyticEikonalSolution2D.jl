export getAnalytic2DeikonalSolution,selfMultiplyWithAnalyticSolution2D,analyticLocal,T02D,getAnalytic2DeikonalSolutionAll

function getAnalytic2DeikonalSolution(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1})

T = zeros(n[1],n[2]);
for j = 1:n[2]
	for i = 1:n[1]
		@inbounds T[i,j] = T02D(h,src,i,j);
	end
end
# finite_vol_integral_for_inverse_r_over_square
Lsrc = 2*(h[1]*asinh(h[2]/h[1]) + h[2]*asinh(h[1]/h[2]))/(h[1]*h[2]);
return T,Lsrc;
end


function T02D(h::Array{Float64,1},src::Array{Int64,1},loc1::Int64,loc2::Int64)
return sqrt( ((loc1 - src[1])*h[1])^2 + ((loc2 - src[2])*h[2])^2 );
end

function selfMultiplyWithAnalyticSolution2D(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1},T::Union{Array{Float32,2},Array{Float64,2}})
for j = 1:n[2]
	for i = 1:n[1]
		@inbounds T[i,j] *= T02D(h,src,i,j);
	end
end
return;
end



function selfMultiplyWithAnalyticSolution2D(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1},T::Union{Array{Float32,1},Array{Float64,1}})
loc = zeros(Int64,2);
for k=1:length(T)
	cs2loc(loc,k,n);
	@inbounds T[k] *= T02D(h,src,loc[1],loc[2]);
end
return;
end




function analyticLocal(h::Array{Float64,1},src::Array{Int64,1},loc1::Int64,loc2::Int64)
@inbounds G01loc = (loc1 - src[1])*h[1];
@inbounds G02loc = (loc2 - src[2])*h[2];
T0loc = sqrt(G01loc*G01loc + G02loc*G02loc);
invt0 = 1/T0loc;
G01loc *= invt0;
G02loc *= invt0;
return T0loc,G01loc,G02loc
end

function analyticLocal(h::Array{Float64,1},src::Array{Int64,1},loc1::Int64,loc2::Int64,GradAns::Array{Float64,1})
@inbounds G01loc = (loc1 - src[1])*h[1];
@inbounds G02loc = (loc2 - src[2])*h[2];
T0loc = sqrt(G01loc*G01loc + G02loc*G02loc);
invt0 = 1/T0loc;
G01loc *= invt0;
G02loc *= invt0;
GradAns[1] = G01loc;
GradAns[2] = G02loc;
return T0loc
end



function getAnalytic2DeikonalSolutionAll(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1})

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];

X1,X2 = ndgrid((0:(n[1]-1))*h[1] - source1,(0:(n[2]-1))*h[2] - source2);

r = sqrt(X1.^2 + X2.^2);
T = r;
L = 1./r;
#Check that for h1=h2:
# finite_vol_integral_for_inverse_r_over_square = 3.52;
# L[src[1],src[2]] = (finite_vol_integral_for_inverse_r_over_square)/sqrt(h[1]*h[2]);
L[src[1],src[2]] = getL0AtSrc2D(h);


G2 = X2.*L;
G2[src[1],src[2]] = 1/sqrt(2);

G1 = X1.*L;
G1[src[1],src[2]] = 1/sqrt(2);

return T,G1,G2,L;
end

function getL0AtSrc2D(h)
	return (2*(h[1]*asinh(h[2]/h[1]) + h[2]*asinh(h[1]/h[2]))/(h[1]*h[2]));
end


## A = 0.25*sqrt(2/Fsrc)*exp(1i*pi/4)*(1/sqrt(pi))*sqrt((1./r));
## const_for_integral = 1.76775; % or 1.76775 for a square cell
## A(src(1),src(2)) = 0.25*sqrt(2/Fsrc)*exp(1i*pi/4)*(1/sqrt(pi))*(const_for_integral/(sqrt(sqrt(h(1)*h(2)))));