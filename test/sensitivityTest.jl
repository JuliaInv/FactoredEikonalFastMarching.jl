using jInv.Mesh;
using FactoredEikonalFastMarching;

include("../examples/getAnalyticalMediums.jl");


function runSensitivityTest()
tic()
println("########### 2D TEST SETTING ##############################");
n = [100,100]
Omega = [0.0,4.0,0.0,8.0];
Msh = getRegularMesh(Omega,n-1);

(kappaSquared,src,T_exact) = getAnalyticalConstGradInv2D(Msh.n+1,Msh.h);
T0 = getAnalytic2DeikonalSolution(n,Msh.h,src)[1];
println("First order test");
doSensTest(kappaSquared,src,T_exact,Msh,T0,false);
println("Second order test");
doSensTest(kappaSquared,src,T_exact,Msh,T0,true);
##########################################################
println("########### 3D TEST SETTING ##############################");
n = [64,64,32]
Omega = [0.0,1.0,0.0,1.0,0.0,0.5];
Msh = getRegularMesh(Omega,n-1);
(kappaSquared,src,T_exact) = getAnalyticalConstGradInv3D(Msh.n+1,Msh.h);
T0 = getAnalytic3DeikonalSolution(n,Msh.h,src)[1];
println("First order test");
doSensTest(kappaSquared,src,T_exact,Msh,T0,false);
println("Second order test");
doSensTest(kappaSquared,src,T_exact,Msh,T0,true);

##########################################################
toc()
return;
end


function doSensTest(kappaSquared,src,T_exact,Msh,T0,HO)
n = Msh.n+1;
pMem = getEikonalTempMemory(n);
pEik = getEikonalParam(Msh,kappaSquared,src,HO);
solveFastMarchingUpwindGrad(pEik,pMem);
T = T0[:].*pEik.T1[:];
println("Maximum error norm: ",maximum(abs((T_exact[:] - T))));
println("Sensitivity Check:");
delta_kappa = 0.05*randn(size(kappaSquared))*mean(kappaSquared[:]);
Tlin = zeros(prod(n));
for k = 1:10
	delta_kappa *= 0.5;
	new_kappa = kappaSquared + delta_kappa;
	pEikNew = getEikonalParam(Msh,new_kappa,src,HO);
	solveFastMarchingUpwindGrad(pEikNew,pMem);
	getSensMatVecEik(delta_kappa[:],Tlin,pEik);
	Tlin += T;
	print("||dm||: ");
	print(norm(delta_kappa[:]))
	print(". ||T(m + dm) - T(m)|| ");
	Tnew = T0[:].*pEikNew.T1[:];
	err0 = (T - Tnew);
	print(norm(err0));
	print(". ||T(m + dm) - (T(m) + J*dm)||: ");
	err = Tnew - Tlin;
	println(norm(err));
	
end
v = rand(prod(n));
Jv = zeros(prod(n));
getSensMatVecEik(v,Jv,pEik);
Z = rand(size(Jv));
JtZ = rand(prod(n));
getSensTMatVecEik(Z,JtZ,pEik);
if abs(dot(JtZ[:],v) - dot(Jv[:],Z[:])) > 1e-8
	println("Checking that (J^T)^T = J for eikonal solver falied: ",abs(dot(JtZ[:],v) - dot(Jv[:],Z[:])) );
end
return;
end
