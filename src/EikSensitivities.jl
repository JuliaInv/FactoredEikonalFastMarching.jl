export getSensMatVecEik, getSensTMatVecEik

function getOp(opCode::Int8,dim::Int64)
if dim==1
	op = mod(opCode,5)-2;
elseif dim==2
	op = mod(opCode,25);
	op = div(op,5)-2;
elseif dim==3
	op = div(opCode,25)-2;
end
end

function getOpCode(op1,op2,op3)
opCode = (op1+2) + 5*(op2+2) + 25*(op3+2);
return opCode;
end

function getSensMatVecEik(b::Array{Float64,1},res::Array{Float64,1},pEik::EikonalParam)
# This is basically a lower triangular solution.
	order = pEik.ordering;
	n = pEik.Mesh.n+1;
	if pEik.Mesh.dim==2
		rowLength = 5;
		if pEik.HO 
			getRowOfSensitivity = getRowOfSensitivity2D;
		else
			getRowOfSensitivity = getRowOfSensitivity2DFO;
		end
		srcCS = loc2cs(pEik.src[1],pEik.src[2],n);
		offset = [1;n[1]];
	else
		rowLength = 7;
		getRowOfSensitivity = getRowOfSensitivity3D;
		srcCS = loc2cs3D(pEik.src[1],pEik.src[2],pEik.src[3],n);
		offset = [1;n[1];n[1]*n[2]];
	end
	G0 = zeros(pEik.Mesh.dim);
	Idxs = zeros(Int64,pEik.Mesh.dim);
	colIdxs = zeros(Int64,rowLength);
	vals 	= zeros(Float64,rowLength);
	############ Invert a lower triangular Matrix ##########################
	res[:] = 0.0;
	for j=order
		r = b[j];
		rowLength_j = getRowOfSensitivity(pEik,j,vals,colIdxs,srcCS,offset,G0,Idxs,n);
		for i = 1:rowLength_j
			r -= vals[i]*res[colIdxs[i]];
		end
		res[j] = r/vals[1];
	end
	if pEik.Mesh.dim==2
		selfMultiplyWithAnalyticSolution2D(n, pEik.Mesh.h,pEik.src,res);
	else
		selfMultiplyWithAnalyticSolution3D(n, pEik.Mesh.h,pEik.src,res);
	end
	
	return res;
end


function getSensTMatVecEik(b::Array{Float64,1},res::Array{Float64,1},pEik::EikonalParam)
# This is basically an upper triangular solution. This function ruins b.
	order = pEik.ordering;
	n = pEik.Mesh.n+1;	
	if pEik.Mesh.dim==2
		rowLength = 5;
		if pEik.HO 
			getRowOfSensitivity = getRowOfSensitivity2D;
		else
			getRowOfSensitivity = getRowOfSensitivity2DFO;
		end
		srcCS = loc2cs(pEik.src[1],pEik.src[2],n);
		offset = [1;n[1]];
	else
		rowLength = 7;
		getRowOfSensitivity = getRowOfSensitivity3D;
		srcCS = loc2cs3D(pEik.src[1],pEik.src[2],pEik.src[3],n);
		offset = [1;n[1];n[1]*n[2]];
	end
	G0 = zeros(pEik.Mesh.dim);
	Idxs = zeros(Int64,pEik.Mesh.dim);
	colIdxs = zeros(Int64,rowLength);
	vals 	= zeros(Float64,rowLength);
	
	############ Invert an upper triangular Matrix ##########################
	N = length(b)
	res[:] = b[:]; # up to index i (backwards, starting from n) res holds the result. From that to 1 it holds the residual.
	if pEik.Mesh.dim==2
		selfMultiplyWithAnalyticSolution2D(n, pEik.Mesh.h,pEik.src,res);
	else
		selfMultiplyWithAnalyticSolution3D(n, pEik.Mesh.h,pEik.src,res);
	end
	
	
	for k=0:(N-1)
		j = order[N-k];
		rowLength_j = getRowOfSensitivity(pEik,j,vals,colIdxs,srcCS,offset,G0,Idxs,n);
		d = res[colIdxs[1]]/vals[1];
		res[colIdxs[1]] = d;
		for i = 2:rowLength_j
			res[colIdxs[i]] -= vals[i]*d;
		end
	end
	return res;
end

function getRowOfSensitivity3D(pEik::EikonalParam,rowIdx::EikIdxType,vals::Array{Float64},colIdxs::Array{Int64},srcCS::Int64,offsetVec::Array{Int64,1},
								G0::Array{Float64,1},Idxs::Array{Int64,1},n::Array{Int64,1});
# the ordering is : [diag,1st cell dx,2nd cell dx,1st cell dy,2nd cell dy]
# in 2D the length of the rows should be of length 5 at most.

OP = pEik.OP;

h = pEik.Mesh.h;
src = pEik.src;
rowIdx = convert(Int64,rowIdx);
colIdxs[1] = rowIdx;
T1 = pEik.T1;
if srcCS==rowIdx	
	vals[1] = 2.0*T1[rowIdx];
	return 1;
end

vals[1] = 0.0;
cs2loc3D(Idxs,rowIdx,n);
@inbounds T0loc = analyticLocal3D(h,src,Idxs[1],Idxs[2],Idxs[3],G0);

k = 1;

for dim = 1:3
	hinvij  = T0loc./h[dim];
	hinv2ij = 0.5*hinvij;
	op = getOp(OP[rowIdx],dim);
	offset = offsetVec[dim];
	if op == 1 # forward 1
		diagDhat  = -hinvij + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinvij*T1[rowIdx+offset]);
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = hinvij*diagSensValue;
		colIdxs[k]  = rowIdx + offset;
	elseif op == -1 # backward 1
		diagDhat  = hinvij + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] - hinvij*T1[rowIdx - offset]);
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = -hinvij*diagSensValue;
		colIdxs[k]  = rowIdx - offset;
	elseif op == -2 # backward 2
		diagDhat  = hinv2ij*3.0 + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinv2ij*(-4.0*T1[rowIdx-offset] + T1[rowIdx-2*offset]));
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = -hinv2ij*4.0*diagSensValue;
		colIdxs[k]  = rowIdx - offset;
		k+=1;
		vals[k]    = hinv2ij*diagSensValue;
		colIdxs[k] = rowIdx - 2*offset;
	elseif op == 2 # forward 2
		diagDhat  = -hinv2ij*3.0 + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinv2ij*(4.0*T1[rowIdx+offset] - T1[rowIdx+2*offset]));
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = hinv2ij*4.0*diagSensValue;
		colIdxs[k]  = rowIdx + offset;
		k+=1;
		vals[k]    = -hinv2ij*diagSensValue;
		colIdxs[k] = rowIdx + 2*offset;
	end
end
return k;
end


function getRowOfSensitivity2D(pEik::EikonalParam,rowIdx::EikIdxType,vals::Array{Float64},colIdxs::Array{Int64},srcCS::Int64,offsetVec::Array{Int64,1},
								G0::Array{Float64,1},Idxs::Array{Int64,1},n::Array{Int64,1});
# the ordering is : [diag,1st cell dx,2nd cell dx,1st cell dy,2nd cell dy]
# in 2D the length of the rows should be of length 5 at most.

OP = pEik.OP;

h = pEik.Mesh.h;
src = pEik.src;

rowIdx = convert(Int64,rowIdx);
colIdxs[1] = rowIdx;
T1 = pEik.T1;
if srcCS==rowIdx	
	vals[1]    = 2.0*T1[rowIdx];
	return 1;
end
vals[1] = 0.0;

cs2loc(Idxs,rowIdx,n);
T0loc = analyticLocal(h,src,Idxs[1],Idxs[2],G0);

k = 1;

for dim = 1:2
	hinvij  = T0loc./h[dim];
	hinv2ij = 0.5*hinvij;
	op = getOp(OP[rowIdx],dim);
	offset = offsetVec[dim];
	if op == 1 # forward 1
		diagDhat  = -hinvij + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinvij*T1[rowIdx+offset]);
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = hinvij*diagSensValue;
		colIdxs[k]  = rowIdx + offset;
	elseif op == -1 # backward 1
		diagDhat  = hinvij + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] - hinvij*T1[rowIdx - offset]);
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = -hinvij*diagSensValue;
		colIdxs[k]  = rowIdx - offset;
	elseif op == -2 # backward 2
		diagDhat  = hinv2ij*3.0 + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinv2ij*(-4.0*T1[rowIdx-offset] + T1[rowIdx-2*offset]));
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = -hinv2ij*4.0*diagSensValue;
		colIdxs[k]  = rowIdx - offset;
		k+=1;
		vals[k]    = hinv2ij*diagSensValue;
		colIdxs[k] = rowIdx - 2*offset;
	elseif op == 2 # forward 2
		diagDhat  = -hinv2ij*3.0 + G0[dim];
		diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinv2ij*(4.0*T1[rowIdx+offset] - T1[rowIdx+2*offset]));
		vals[1]                 += diagDhat*diagSensValue;
		k+=1;
		vals[k]     = hinv2ij*4.0*diagSensValue;
		colIdxs[k]  = rowIdx + offset;
		k+=1;
		vals[k]    = -hinv2ij*diagSensValue;
		colIdxs[k] = rowIdx + 2*offset;
	end
end
return k;
end


function getRowOfSensitivity2DFO(pEik::EikonalParam,rowIdx::EikIdxType,vals::Array{Float64},colIdxs::Array{Int64},srcCS::Int64,offsetVec::Array{Int64,1},
								G0::Array{Float64,1},Idxs::Array{Int64,1},n::Array{Int64,1});
# the ordering is : [diag,1st cell dx,2nd cell dx,1st cell dy,2nd cell dy]
# in 2D the length of the rows should be of length 5 at most.

OP  = pEik.OP;
h   = pEik.Mesh.h;
src = pEik.src;
T1  = pEik.T1;
rowIdx = convert(Int64,rowIdx);

@inbounds colIdxs[1] = rowIdx;

if srcCS==rowIdx	
	@inbounds vals[1]    = 2.0*T1[rowIdx];
	return 1;
end
@inbounds vals[1] = 0.0;

cs2loc(Idxs,rowIdx,n);

@inbounds T0loc = analyticLocal(h,src,Idxs[1],Idxs[2],G0);

k = 1;

for dim = 1:2
	@inbounds op = getOp(OP[rowIdx],dim);
	if op != 0
		@inbounds firstNeighbor = rowIdx + op*offsetVec[dim];
		@inbounds hinvij  = op*(T0loc./h[dim]);
	
		@inbounds diagDhat  = G0[dim]-hinvij;
		@inbounds diagSensValue = 2.0*(diagDhat*T1[rowIdx] + hinvij*T1[firstNeighbor]);
		@inbounds vals[1]     += diagDhat*diagSensValue;
		k+=1;
		@inbounds vals[k]     = hinvij*diagSensValue;
		@inbounds colIdxs[k]  = firstNeighbor;
	end
end
return k;
end
