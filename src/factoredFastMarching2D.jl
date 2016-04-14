function factoredFastMarching(pEik::EikonalParam, mem::EikonalTempMemory)
# time1 = 0.0;
# time2 = 0.0;
# time3 = 0.0;


kappaSquared 	= pEik.kappaSquared;
h 				= pEik.Mesh.h;
src 			= pEik.src;
n 				= pEik.Mesh.n+1;
HO 				= pEik.HO;

if isempty(pEik.T1)
	pEik.T1 = zeros(Float64,tuple(n...));
end
T 		 = pEik.T1;
if isempty(pEik.ordering)
	pEik.ordering = zeros(EikIdxType,prod(n));
end
ordering = pEik.ordering;
if isempty(pEik.OP)
	pEik.OP = zeros(Int8,tuple(n...));
end
OP       = pEik.OP;
OP[:] 	 = 0;



T[:] = Inf;
Done = mem.Done;
Done[:] = false;

heapArr = mem.V;
heapIdx = mem.J;
frontHeap = initHeap(heapArr,heapIdx);

monoThresh = max(h[1],h[2])^2;


# Done = zeros(Bool,n[1]+4,n[2]+4); # the +4 is actually only for HO. For FO 2 is enough but need to change access to Done in helpFuncs




src_cs = loc2cs(src,n);

kappaSrc = sqrt(kappaSquared[src_cs]);
T[src_cs] = kappaSrc;
insertToHeap(frontHeap,0.0,src_cs);

k_order = 1;
Neighbors2D = zeros(Int64,4,2);
Boundary2D = zeros(Bool,4);
curr_loc = zeros(Int64,2);
GradLoc = zeros(Float64,2);
A = zeros(Float64,2);
B = zeros(Float64,2);
OpLoc = zeros(Int8,2);
curr_val = zeros(Float64,1);

hinv = 1./h;

while frontHeap.size > 0
	# tic();
	curr = getMin(frontHeap,curr_val);
	cs2loc(curr_loc,curr,n);
	while getPaddedArr2D(Done,curr_loc[1],curr_loc[2]) && frontHeap.size>0
		curr = getMin(frontHeap,curr_val);
		cs2loc(curr_loc,curr,n);
	end
	# time1 += toq();
	if getPaddedArr2D(Done,curr_loc[1],curr_loc[2]) && frontHeap.size==0
        break;
    end
	ordering[k_order] = curr;
    k_order += 1;
	setPaddedArr2D(Done,curr_loc[1],curr_loc[2],true); 
	get2DNeighbors(curr_loc,n,Neighbors2D,Boundary2D);
	for k = 1:4
		if Boundary2D[k]
			@inbounds nei_loc1 = Neighbors2D[k,1];
			@inbounds nei_loc2 = Neighbors2D[k,2];
			
			if !getPaddedArr2D(Done,nei_loc1,nei_loc2) 
				# tic();
				T0loc = analyticLocal(h,src,nei_loc1,nei_loc2,GradLoc);
				
				T_k = getUpdatedFactoredVal_2D(HO,kappaSquared,T,Done,nei_loc1,nei_loc2,h,hinv,src,T0loc,GradLoc,OpLoc,A,B);
				
				# @inbounds if T_k*T0loc < curr_val[1]
					# @inbounds if abs(T_k*T0loc-curr_val[1]) > monoThresh
						 # println("Monotonicity not fulfilled by ",(T_k*T0loc-curr_val)," - not corrected. HO = ",HO);
					# end
				# end
				# time2 += toq();
				# tic();
				
				if T_k < T[nei_loc1,nei_loc2]
					@inbounds T[nei_loc1,nei_loc2] = T_k;
					@inbounds OP[nei_loc1,nei_loc2] = getOpCode(OpLoc[1],OpLoc[2],0) ;
					@inbounds insertToHeap(frontHeap,T0loc*T_k,loc2cs(nei_loc1,nei_loc2,n));
				end
				# time3 += toq();
			end
			
		end
	end
end

# println("time heap: ",time1, " time quad: ",time2, " time heapInsert: ",time3);
return T,ordering;
end;