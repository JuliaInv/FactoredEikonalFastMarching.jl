
export getAnalytic3DeikonalSolution,selfMultiplyWithAnalyticSolution3D,analyticLocal3D,T03D,getAnalytic3DeikonalSolutionAll


function getAnalytic3DeikonalSolution(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1})
T = zeros(n[1],n[2],n[3]);

for k = 1:n[3]
	for j = 1:n[2]
		for i = 1:n[1]
			T[i,j,k] = T03D(n,h,src,i,j,k);
		end
	end
end
return T,0.0;
end




function selfMultiplyWithAnalyticSolution3D(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1},T::Union{Array{Float32,3},Array{Float64,3}})
for k = 1:n[3]
	for j = 1:n[2]
		for i = 1:n[1]
			T[i,j,k] *= T03D(n,h,src,i,j,k);
		end
	end
end
return;
end



function selfMultiplyWithAnalyticSolution3D(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1},T::Union{Array{Float32,1},Array{Float64,1}})
loc = zeros(Int64,3);
for k=1:length(T)
	cs2loc3D(loc,k,n);
	T[k] *= T03D(n,h,src,loc[1],loc[2],loc[3]);
end
return;
end


function T03D(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1},loc1::Int64,loc2::Int64,loc3::Int64)
return sqrt(((loc1 - src[1])*h[1])^2 + ((loc2 - src[2])*h[2])^2 + ((loc3 - src[3])*h[3])^2);
end

function analyticLocal3D(h::Array{Float64,1},src::Array{Int64,1},loc1::Int64,loc2::Int64,loc3::Int64)
G01loc = (loc1 - src[1])*h[1];
G02loc = (loc2 - src[2])*h[2];
G03loc = (loc3 - src[3])*h[3];

T0loc = sqrt(G01loc*G01loc + G02loc*G02loc + G03loc*G03loc);
invt0 = 1/T0loc;
G01loc *= invt0;
G02loc *= invt0;
G03loc *= invt0;
return T0loc,G01loc,G02loc,G03loc
end

function analyticLocal3D(h::Array{Float64,1},src::Array{Int64,1},loc1::Int64,loc2::Int64,loc3::Int64,GradAns::Array{Float64,1})
G01loc = (loc1 - src[1])*h[1];
G02loc = (loc2 - src[2])*h[2];
G03loc = (loc3 - src[3])*h[3];

T0loc = sqrt(G01loc*G01loc + G02loc*G02loc + G03loc*G03loc);
invt0 = 1/T0loc;
G01loc *= invt0;
G02loc *= invt0;
G03loc *= invt0;
GradAns[1] = G01loc;
GradAns[2] = G02loc;
GradAns[3] = G03loc;
return T0loc
end


function getAnalytic3DeikonalSolutionAll(n::Array{Int64,1},h::Array{Float64,1},src::Array{Int64,1})

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];
source3 = (src[3]-1)*h[3];

X1,X2,X3 = ndgrid((0:(n[1]-1))*h[1] - source1,(0:(n[2]-1))*h[2] - source2,(0:(n[3]-1))*h[3] - source3);

r = sqrt(X1.^2 + X2.^2 + X3.^2);
T = r;
L = 1./r;
G2 = X2.*L;
G2[src[1],src[2],src[3]] = 1/sqrt(3);

G1 = X1.*L;
G1[src[1],src[2],src[3]] = 1/sqrt(3);

G3 = X3.*L;
G3[src[1],src[2],src[3]] = 1/sqrt(3);

L[src[1],src[2],src[3]] = getL0AtSrc3D(h);

return T,G1,G2,G3,L;
end


function getL0AtSrc3D(h)
# finite_vol_integral_for_inverse_r_over_cube
# Solving the 3D integral for the source of the laplacian:

# L(src) = (1/vol)*int_{-h1/2<x1<h1/2,-h2/2<x2<h2/2,-h3/2<x3<h3/2}{1/sqrt(x^2 + y^2 + z^2)}  
# vol = h1*h2*h3;
# Now we reduce one dimension from the integration and use 2D numerical mid-point integration.
# from quickmath: int_{-h/2<x<h/2}(1/sqrt(x^2 + y^2)) = 2*asinh(h/(2*abs(y)));

# Therefore: we reduce the first coordinate since we expect the third one to be the smallest

###########################################################################
# tic()
# n_num = 100;
# x1_num = -h[1]/2+h[1]/(2*n_num):h[1]/n_num:h[1]/2-h[1]/(2*n_num);
# x2_num = -h[2]/2+h[2]/(2*n_num):h[2]/n_num:h[2]/2-h[2]/(2*n_num);
# x3_num = -h[3]/2+h[3]/(2*n_num):h[3]/n_num:h[3]/2-h[3]/(2*n_num);

# X1sq,X2sq,X3sq = ndgrid(x1_num.^2,x2_num.^2,x3_num.^2);
# F = 0.0;
# for ii = 1:n_num
	# for jj = 1:n_num
		# for kk = 1:n_num
			# F += 1/sqrt(X1sq[ii,jj,kk]+X2sq[ii,jj,kk]+X3sq[ii,jj,kk]);
		# end
	# end
# end
# size_numerical_box = (h[1]*h[2]*h[3])/(n_num*n_num*n_num)
# F *= size_numerical_box;

# F *= (1/prod(h)); # finite volume definition for the cell.
# println(F)
# toc()
###################################################################
# n_num = 200;
# x2_num = -h[2]/2+h[2]/(2*n_num):h[2]/n_num:h[2]/2-h[2]/(2*n_num);
# x3_num = -h[3]/2+h[3]/(2*n_num):h[3]/n_num:h[3]/2-h[3]/(2*n_num);

# X2sq,X3sq = ndgrid(x2_num.^2,x3_num.^2);
# F = 0.0;
# for ii = 1:n_num
	# for jj = 1:n_num
		# F += 2*asinh((h[1])./(2*sqrt(X2sq[ii,jj] + X3sq[ii,jj])));
	# end
# end
# size_numerical_box = (h[2]*h[3])/(n_num*n_num)
# F *= size_numerical_box;
# F *= (1/prod(h)); # finite volume definition for the cell.
# println(F);

n_num = 200;
x2_num = -h[2]/2+h[2]/(2*n_num):h[2]/n_num:h[2]/2-h[2]/(2*n_num);
x3_num = -h[3]/2+h[3]/(2*n_num):h[3]/n_num:h[3]/2-h[3]/(2*n_num);

X2sq,X3sq = ndgrid(x2_num.^2,x3_num.^2);
F = 2.0*sum(asinh((0.5*h[1])./sqrt(X2sq+X3sq)));
# F = 0.0;
# for ii = 1:n_num
	# for jj = 1:n_num
		# F += 2*asinh((h[1])./(2*sqrt(X2sq[ii,jj] + X3sq[ii,jj])));
	# end
# end
size_numerical_box = (h[2]*h[3])/(n_num*n_num)
F *= size_numerical_box;
F *= 3*(1/prod(h)); # finite volume definition for the cell.
# println(F)

# F =  (2*h[1]*h[2] + 2*h[2]*h[3] + 2*h[1]*h[3])/(h[1]*h[2]*h[3]);
return F;

end


