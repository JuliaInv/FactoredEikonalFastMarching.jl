function runExperimentAndWriteResults(kappaSquared::Array{Float64},h::Array{Float64,1},src::Array{Int64,1},n::Array{Int64,1},T_exact::Array{Float64},WU::Float64 )
Omega = zeros(2*length(n));
Omega[2:2:end] = (n-1).*h;
Mesh = getRegularMesh(Omega,n-1);
##########################################################

# if length(n)==2
	# matshow(sqrt(kappaSquared)); colorbar();
	# xlabel("y");
	# ylabel("x");
	# xticks(0:div(n[2],8):n[2],0:1:8);
	# yticks(0:div(n[1],4):n[1],0:1:4);
# end


pMem = getEikonalTempMemory(n);
pEik = getEikonalParam(Mesh,kappaSquared,src,false);
pEik.T1 = zeros(Float64,tuple(n...));

tictimes = 1000000;
times = div(tictimes,prod(n)) + 1;
println("Timing is based on averaged ",times," runs.");
tic()
for k=1:times
	solveFastMarchingUpwindGrad(pEik,pMem);
end
t1 = toq();
t1 = t1/times;


if Mesh.dim==2
	selfMultiplyWithAnalyticSolution2D(Mesh.n+1,Mesh.h,src,pEik.T1);
else
	selfMultiplyWithAnalyticSolution3D(Mesh.n+1,Mesh.h,src,pEik.T1)
end



T1 = copy(pEik.T1)

pEik.HO = true;
ERR1 = T_exact - pEik.T1;

tic()
for k=1:times	
	solveFastMarchingUpwindGrad(pEik,pMem);
end
t2 = toq()
t2 = t2/times;


if Mesh.dim==2
	selfMultiplyWithAnalyticSolution2D(Mesh.n+1,Mesh.h,src,pEik.T1);
else
	selfMultiplyWithAnalyticSolution3D(Mesh.n+1,Mesh.h,src,pEik.T1)
end


sqrtN = sqrt(prod(n));
T2 = copy(pEik.T1)

ERR2 = T_exact - pEik.T1;

@printf("[%3.4f,%3.4f]\t[%d,%d]\t[%3.2e,%3.2e]\t%3.2fs(%3.2f)\t",h[1],h[2],n[1],n[2],vecnorm(ERR1,Inf),vecnorm(ERR1,2)/sqrtN,t1,t1/WU);

@printf("[%3.2e,%3.2e]\t%3.2fs(%3.2f)\n",vecnorm(ERR2,Inf),vecnorm(ERR2,2)/sqrtN,t2,t2/WU);

# if length(n)==2
	# figure()
	# contour(T_exact[end:-1:1,:],20,linestyles = "-",colors = "black",linewidths=1.5)
	# contour(T1[end:-1:1,:],20,linestyles = ":",colors = "red", linewidths=1.5)
	# contour(T2[end:-1:1,:],20,linestyles = "--",colors = "blue", linewidths=1.5)
	# xlabel("y");
	# ylabel("x");
	# xticks(0:div(n[2],8):n[2],0:1:8);
	# yticks(0:div(n[1],4):n[1],4:-1:0);

	
	# matshow(abs(T1 - T_exact)); colorbar();
	# xlabel("y");
	# ylabel("x");
	# xticks(0:div(n[2],8):n[2],0:1:8);
	# yticks(0:div(n[1],4):n[1],0:1:4);
	
	# matshow(abs(T2 - T_exact)); colorbar();
	# xlabel("y");
	# ylabel("x");
	# xticks(0:div(n[2],8):n[2],0:1:8);
	# yticks(0:div(n[1],4):n[1],0:1:4);
# end

# if length(n)==2
	# figure()
	# contour(T_exact[end:-1:1,:],100,linestyles = "-",colors = "black",linewidths=2.0)
	# contour(T1[end:-1:1,:],100,linestyles = ":",colors = "red", linewidths=2.0)
	# contour(T2[end:-1:1,:],100,linestyles = "--",colors = "blue", linewidths=2.0)
	# xlabel("y");
	# ylabel("x");
	# xticks(0:div(n[2],16):n[2],0:0.5:8);
	# yticks(0:div(n[1],40):n[1],4:-0.1:0);

	
# end


return;
end