type minHeap
    valArr::Array{Float64,1};
	indArr::Array{Int64,1};
    size::Int64;
    maxSize::Int64;
end

function initHeap(n::Int64)
valArr = zeros(n);
indArr = zeros(Int64,n);
heap = minHeap(valArr,indArr,0,n);
return heap;
end

function initHeap(valArr::Array{Float64,1},indArr::Array{Int64,1})
if length(valArr)!=length(indArr)
	error("initHeap::unevenly lengthed arrays");
end
valArr[:] = 0.0;
indArr[:] = 0;
heap = minHeap(valArr,indArr,0,length(valArr));
return heap;
end

function insertToHeap(heap::minHeap,val::Float64,ind::Int64)
if heap.size == heap.maxSize
	println("Allocating more memory in heap.")
	tmpArr = heap.valArr;
	heap.valArr = zeros(2*heap.maxSize);
	heap.valArr[1:heap.maxSize] = tmpArr;
	tmpArr = heap.indArr;
	heap.indArr = zeros(Int64,2*heap.maxSize);
	heap.indArr[1:heap.maxSize] = tmpArr;
	tmpArr = 0; # releasing memory???
	heap.maxSize = 2*heap.maxSize;
end
heap.size += 1;
now = heap.size;
heap.valArr[now] = val;
heap.indArr[now] = ind;


# par = div(now,2);
# while now > 1 && val < heap.valArr[par]
    # # @inbounds if ()
    # @inbounds heap.valArr[now] = heap.valArr[par];
	# @inbounds heap.indArr[now] = heap.indArr[par];
    # @inbounds heap.valArr[par] = val;
	# @inbounds heap.indArr[par] = ind;
    # now = par;
	# par = div(now,2);
    # # else
        # # break;
    # # end
# end




while now > 1
    par = div(now,2);
    @inbounds if (val < heap.valArr[par])
        @inbounds heap.valArr[now] = heap.valArr[par];
		@inbounds heap.indArr[now] = heap.indArr[par];
        @inbounds heap.valArr[par] = val;
		@inbounds heap.indArr[par] = ind;
        now = par;
    else
        break;
    end
end
return;
end


function getMin(heap::minHeap)
if heap.size >= 1
	if (heap.size==0)
		least = NaN;
		ind = 0;
		return least,ind;
	elseif heap.size==1
		least = heap.valArr[1];
		ind = heap.indArr[1];
		heap.valArr[1] = 0.0;
		heap.indArr[1] = 0;
		heap.size = 0;
		return least,ind;
	end
end

# taking the root node as the minimum
@inbounds least = heap.valArr[1];
@inbounds least_ind = heap.indArr[1];
# setting the last element as x and reducing the size of the heap.
@inbounds x = heap.valArr[heap.size];
@inbounds x_ind = heap.indArr[heap.size];
@inbounds heap.valArr[heap.size] = 0.0;
@inbounds heap.indArr[heap.size] = 0;
heap.size = heap.size - 1;

curr = 1;

while true
    leftChild = 2*curr; # 2k is the left child of k 
    if leftChild > heap.size # no left child.
		break;
	else
		rightChild = curr*2 + 1;
        @inbounds if rightChild > heap.size || heap.valArr[leftChild] < heap.valArr[rightChild]
            child = leftChild;
        else
            child = rightChild;
        end
        if x > heap.valArr[child]
            @inbounds heap.valArr[curr] = heap.valArr[child];
		    @inbounds heap.indArr[curr] = heap.indArr[child];
            curr = child;
        else
            break;
        end
    end
end

@inbounds heap.valArr[curr] = x;
@inbounds heap.indArr[curr] = x_ind;
return least,least_ind;
end



function getMin(heap::minHeap,Ans::Array{Float64,1})
if heap.size >= 1
	if (heap.size==0)
		least = NaN;
		@inbounds Ans[1] = least;
		return 0;
	elseif heap.size==1
		Ans[1] = heap.valArr[1];
		ind = heap.indArr[1];
		heap.valArr[1] = 0.0;
		heap.indArr[1] = 0;
		heap.size = 0;
		return ind;
	end
end

# taking the root node as the minimum
@inbounds least = heap.valArr[1];
@inbounds least_ind = heap.indArr[1];
# setting the last element as x and reducing the size of the heap.
@inbounds x = heap.valArr[heap.size];
@inbounds x_ind = heap.indArr[heap.size];
@inbounds heap.valArr[heap.size] = 0.0;
@inbounds heap.indArr[heap.size] = 0;
heap.size = heap.size - 1;

curr = 1;
# leftChild = 2; 
# while leftChild <= heap.size # no left child. 
   # rightChild = curr*2 + 1;
   # child = rightChild + (heap.valArr[leftChild] < heap.valArr[rightChild] || rightChild > heap.size)*(leftChild - rightChild)
   
   # # @inbounds if heap.valArr[leftChild] < heap.valArr[rightChild] || rightChild > heap.size
      # # child = leftChild;
   # # else
      # # child = rightChild;
   # # end
   # if x > heap.valArr[child]
       # @inbounds heap.valArr[curr] = heap.valArr[child];
	   # @inbounds heap.indArr[curr] = heap.indArr[child];
       # curr = child;
   # else
       # break;
   # end
   # leftChild = 2*curr; # 2k is the left child of k 
# end


while true
    leftChild = 2*curr; # 2k is the left child of k 
    if leftChild > heap.size # no left child.
		break;
	else
		rightChild = curr*2 + 1;
        @inbounds if rightChild > heap.size || heap.valArr[leftChild] < heap.valArr[rightChild]
            child = leftChild;
        else
            child = rightChild;
        end
        if x > heap.valArr[child]
            @inbounds heap.valArr[curr] = heap.valArr[child];
		    @inbounds heap.indArr[curr] = heap.indArr[child];
            curr = child;
        else
            break;
        end
    end
end

@inbounds heap.valArr[curr] = x;
@inbounds heap.indArr[curr] = x_ind;
@inbounds Ans[1] = least;
return least_ind;
end

