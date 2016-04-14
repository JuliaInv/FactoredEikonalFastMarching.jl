function factoredFastMarching3D(pEik::EikonalParam, mem::EikonalTempMemory)

kappaSquared 	= pEik.kappaSquared;
h 				= pEik.Mesh.h;
src 			= pEik.src;
n 				= pEik.Mesh.n+1;
HO 				= pEik.HO;

N = prod(n);

if isempty(pEik.T1)
	pEik.T1 = zeros(Float64,tuple(n...));
end
T = pEik.T1;
if isempty(pEik.ordering)
	pEik.ordering = zeros(EikIdxType,N);
end

ordering = pEik.ordering;
if isempty(pEik.OP)
	pEik.OP = zeros(Int8,tuple(n...));
end
OP       = pEik.OP;
OP[:] 	 = 0;


Done = mem.Done;
Done[:] = false;

heapArr = mem.V;
heapIdx = mem.J;
frontHeap = initHeap(heapArr,heapIdx);

# time1 = 0.0;
# time2 = 0.0;
# time3 = 0.0;

monoThresh = max(h[1],h[2],h[3])^2;
T[:] = Inf;

hinv = 1./h;
src_cs = loc2cs3D(src,n);

kappaSrc = sqrt(kappaSquared[src_cs]);
T[src_cs] = kappaSrc;

insertToHeap(frontHeap,0.0,src_cs);
k_order = 1;
Neighbors3D = zeros(Int64,6,3);
notBoundary3D = zeros(Bool,6);
curr_loc = zeros(Int64,3);
GradLoc = zeros(Float64,3);
A = zeros(Float64,3);
B = zeros(Float64,3);
OpLoc = zeros(Int8,3);
curr_val = zeros(Float64,1);

while frontHeap.size > 0
	# tic();
	curr = getMin(frontHeap,curr_val);
	cs2loc3D(curr_loc,curr,n);
	@inbounds while getPaddedArr3D(Done,curr_loc[1],curr_loc[2],curr_loc[3]) && frontHeap.size>0
		curr = getMin(frontHeap,curr_val);
		cs2loc3D(curr_loc,curr,n);
	end
	@inbounds if getPaddedArr3D(Done,curr_loc[1],curr_loc[2],curr_loc[3]) && frontHeap.size==0
        break;
    end
	# time1 += toq();
	ordering[k_order] = curr;
    k_order += 1;
	@inbounds setPaddedArr3D(Done,curr_loc[1],curr_loc[2],curr_loc[3],true);
	@inbounds get3DNeighbors(curr_loc,n,Neighbors3D,notBoundary3D);
	
	for k = 1:6
		@inbounds if notBoundary3D[k]
			@inbounds nei_loc1 = Neighbors3D[k,1];
			@inbounds nei_loc2 = Neighbors3D[k,2];
			@inbounds nei_loc3 = Neighbors3D[k,3];
			if !getPaddedArr3D(Done,nei_loc1,nei_loc2,nei_loc3) 
				# tic();
				T0loc = analyticLocal3D(h,src,nei_loc1,nei_loc2,nei_loc3,GradLoc);
				m_loc = kappaSquared[nei_loc1,nei_loc2,nei_loc3];
				@inbounds T_k = getUpdatedFactoredVal_3D(HO,m_loc,T,Done,nei_loc1,
									nei_loc2,nei_loc3,h,hinv,n,src,T0loc,GradLoc,OpLoc,A,B);
				
				# @inbounds if T_k*T0loc < curr_val[1]
					# @inbounds if abs(T_k*T0loc-curr_val[1]) > monoThresh
						# println("Monotonicity not fulfilled by ",(T_k*T0loc-curr_val),", not corrected.");
					# end
				# end
				# time2 += toq();
				# tic();
				@inbounds if T_k < T[nei_loc1,nei_loc2,nei_loc3]
					@inbounds T[nei_loc1,nei_loc2,nei_loc3] = T_k;
					@inbounds OP[nei_loc1,nei_loc2,nei_loc3] = getOpCode(OpLoc[1],OpLoc[2],OpLoc[3]) ;
					@inbounds insertToHeap(frontHeap,T0loc*T_k,loc2cs3D(nei_loc1,nei_loc2,nei_loc3,n));
				end
				# time3 += toq();
			end
		end
	end
end
# println("time heap: ",time1, " time quad: ",time2, " time heapInsert: ",time3);
return T,ordering;
end;
