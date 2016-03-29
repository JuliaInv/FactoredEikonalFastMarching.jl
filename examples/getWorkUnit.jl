function GetWorkunitForLoop(n::Array{Int64,1},h::Array{Float64,1})
hinv = 1./h;
times = 20;
t_wu = 0.0;
if length(n)==2
	t = ones(n[1]+2,n[2]+2);
	z = zeros(n[1]+2,n[2]+2);
	tic()
	for k=1:times
		for j = 1:n[2]
			for i = 1:n[1]
				i_t = i+1;
				j_t = j+1;
				z[i,j] = ((0.5*hinv[1]*(t[i_t+1,j_t] - t[i_t-1,j_t]))^2 + (0.5*hinv[2]*(t[i_t,j_t+1] - t[i_t,j_t-1]))^2 );
			end
		end
	end
	t_wu = toq();
elseif length(n)==3
	t = ones(n[1]+2,n[2]+2,n[3]+2);
	z = zeros(n[1]+2,n[2]+2,n[3]+2);
	tic()
	for k=1:times
		for k = 1:n[3]
			for j = 1:n[2]
				for i = 1:n[1]
					i_t = i+1;
					j_t = j+1;
					k_t = k+1;
					z[i,j,k] = ((0.5*hinv[1]*(t[i_t+1,j_t,k_t] - t[i_t-1,j_t,k_t]))^2 + (0.5*hinv[2]*(t[i_t,j_t+1,k_t] - t[i_t,j_t-1,k_t]))^2 + (0.5*hinv[3]*(t[i_t,j_t,k_t+1] - t[i_t,j_t,k_t-1]))^2);
				end
			end
		end
	end
	t_wu = toq();
else
	error("unkonwn dimention");
end
WU = t_wu/times;
println("workunit is: ",WU," for n = ",n);
return WU;
end

# function showWorkunit(h::Array{Float64,1},WU::Float64)
# I = [1.6,1.6,0.8];

# n = zeros(Int64,3);
# n[1] = round(Int64,(I[1]/h[1])+1);
# n[2] = round(Int64,(I[2]/h[2])+1);
# n[3] = round(Int64,(I[3]/h[3])+1);

# t = rand(prod(n));
# (D1,D2,D3) = getLongDiffGradOperators3D(h,n);
# times = 30;
# t_wu = 0.0;
# tic()
# for k=1:times
	# y = D1*t;
	# y.*=y;
	# z = D2*t; 
	# z.*=z;
	# y+=z;
	# z = D3*t; 
	# z.*=z;
	# y+=z;
# end
# t_wu = toq();

# WU = t_wu/times;
# D1 = D2 = D3 = 0.0;
# println("workunit is: ",WU," for n = ",n);
# return;
# end