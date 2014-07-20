
# a DisjointSet is a structure that stores counts of the number of observations that
# fall into each element of the power set of the given sets. Note that a PowerSet
# would be similar, but (for example) the slot for A would contain A&B, whereas for
# a DisjointSet, the slot for A is for observations of (A and not-B and not-C), etc.
type DisjointSet
	labels
	counts

	function DisjointSet(dat::AbstractMatrix, labs::AbstractVector)
		# dat must be 2-D and the num of cols must be equal to the length of the labels
		@assert size(dat,2) == size(labs,1)
		nrow = size(dat,1); ncol = size(dat,2)

		pssize = 2^ncol
		cts = zeros(Int64, pssize)

		for row in 1:nrow
			# which cell to increment? start with 0 and add powers of 2 from right to left
			# then add 1 because Julia is 1-based
			i = 0
			for col in 1:ncol
				if (dat[row, col])
					i += 2^(ncol-col)
				end
			end
			cts[i+1] += 1
		end
		new(labs, cts)
	end
end

# an EulerState is just a Vector{Float64}
typealias EulerState Vector{Float64}

type EulerObject
	nparams
	labels
	lb
	ub
	sizes
	target
	specs
	evalfn
end

type EulerSpec
	shape
	clamp
	statepos

	function EulerSpec(shape, clamp, statepos)
		# this constructor enforces invariants
		if in(shape, [:circle, :square, :triangle])
			@assert length(clamp) == 2
			@assert length(statepos) == 2
		elseif in(shape, [:rectangle])
			@assert length(clamp) == 3
			@assert length(statepos) == 3
		else 
			error("Unknown EulerSpec shape: ", string(shape))
		end
		new(shape, clamp, statepos)
	end
end

EulerSpec() = EulerSpec(:circle)
function EulerSpec(shape)
	if in(shape, [:circle, :square, :triangle])
		EulerSpec(shape, [NaN, NaN], [0, 0])
	elseif in(shape, [:rectangle])
		EulerSpec(shape, [NaN, NaN, NaN], [0, 0, 0])
	else
		error("Unknown EulerSpec shape: ", string(shape))
	end
end

