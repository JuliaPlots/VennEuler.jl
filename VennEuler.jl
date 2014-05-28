
# throwaway file to explore things...

randdata = randbool(20, 3); # 3 cols
setlabels = ["A", "B", "C"];

# a DisjointSet is a structure that stores counts of the number of observations that
# fall into each element of the power set of the given sets. Note that a PowerSet
# would be similar, but (for example) the slot for A would contain A&B, whereas for
# a DisjointSet, the slot for A is for observations of A and not-B, not-C, etc.
type DisjointSet
	labels
	counts

	function DisjointSet(dat::AbstractArray, labs::AbstractArray)
		# dat must be 2-D and the num of cols must be equal to the length of the labels
		@assert length(size(dat)) == 2
		@assert length(size(labs)) == 1
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

# pretty print a DisjointSet

# regress two DisjointSets

# convert a DisjointSet to a PowerSet -- needed?

type EulerState
	labels::Vector{ASCIIString}
	elems::Vector{EulerStateElem}
end

abstract EulerStateElem

type EulerStateCircle <: EulerStateElem
	area
	params
	lowerbounds
	upperbounds
end

# method to generate an EulerState from a set of label names

# method to generate a 2-D bitmap from an EulerStateElem

# method to generate all bitmaps from an EulerState

# method to compute overlaps from the bitmaps from an EulerState, yielding a DisjointSet



