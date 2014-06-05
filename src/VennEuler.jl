module VennEuler

using Base
using Iterators
using NLopt
using Cairo

export
	DisjointSet,
	EulerState,
	EulerObject,
	makeeulerobject,
	evaleulerstate,
	optimize,
	render


# a DisjointSet is a structure that stores counts of the number of observations that
# fall into each element of the power set of the given sets. Note that a PowerSet
# would be similar, but (for example) the slot for A would contain A&B, whereas for
# a DisjointSet, the slot for A is for observations of A and not-B, not-C, etc.
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

# pretty print a DisjointSet

# regress two DisjointSets

# convert a DisjointSet to a PowerSet -- needed?


# an EulerState is just a Vector{Float64}
typealias EulerState Vector{Float64}

# to create one, a function takes specifications and returns several things:
# an initial EulerState
# a closure that can be used to generate scores for an EulerState
# lower bound vector, upper bound vector
# areas?

type EulerObject
	nparams
	labels
	lb
	ub
	sizes
	target
	evalfn
end

dupeelements(qq) = [qq[ceil(i)] for i in (1:(2*length(qq)))/2]

function makeeulerobject(labels::Vector{ASCIIString}, sizes::Vector, target::DisjointSet)
	@assert length(labels) == length(sizes)
	# create the bounds vectors
	nparams = 2 * length(labels)
	# normalize the sizes array so that the sum is = 1/2 
	sizes = sizes / (2 * sum(sizes))
	# this forces the centers to not overlap with the boundary
	lb = dupeelements(sizes)
	ub = dupeelements(1 .- sizes)
	@assert all(lb .< ub)
	# create the initial state vector
	es::EulerState = rand(nparams) .* (ub .- lb) .+ lb
	
	# return: state vector, state object (with bounds, closure, etc)
	eo = EulerObject(nparams, labels, lb, ub, sizes, target, identity)
	eo.evalfn = (x,g) -> begin
		@assert length(g) == 0
		cost = VennEuler.evaleulerstate(eo, x) 
		#println(x, " (", cost, ")")
		cost
	end
	(es, eo)
end

function optimize(obj::EulerObject, state::EulerState; xtol=1/200, ftol=1/1000, maxtime=5, init_step=.1)
	opt = Opt(:GN_DIRECT, length(state))
	lower_bounds!(opt, obj.lb)
	upper_bounds!(opt, obj.ub)
	xtol_abs!(opt, xtol)
	ftol_abs!(opt, ftol)
	min_objective!(opt, obj.evalfn)
	maxtime!(opt, maxtime)
	initial_step!(opt, init_step)
	NLopt.optimize(opt, state)
end

function evaleulerstate(obj::EulerObject, state::EulerState; verbose::Int64=0, px::Int64=200)
	# given this state vector and the object, do the following:
	# generate a 2-D bitmap from each object
	# foreach element of the DisjointSet, calculate the size of the overlap of the bitmaps
	# compare the overlaps with the target, returning the error metric
	bitmaps = [makebitmapcircle(state[2i-1], state[2i], obj.sizes[i], px) 
				for i in 1:length(obj.labels)]

	# iterate through the powerset index
	# ignore 000
	overlaps = zeros(2^length(obj.labels)-1)
	for psi = 1:length(overlaps)
		bpsi = bits(psi)
		compare_str = bpsi[(endof(bpsi) - length(obj.labels) + 1):end]
		if (verbose > 0)
			@show compare_str
		end
		# if we have 101, then take the overlap of bitmaps[1], !bitmaps[2], and bitmaps[3]
		overlap = trues(px,px) 
		for i in 1:length(compare_str)
			overlap = overlap & (compare_str[i] == '1' ? bitmaps[i] : !bitmaps[i])
		end
		if (verbose > 0) showbitmap(overlap) end
		overlaps[psi] = sum(overlap)
	end
	if (verbose > 0)
		@show overlaps
	end
	# TODO: use better error metric
	# compare overlaps with obj.target
	overlaps_norm = overlaps ./ sum(overlaps)
	target_norm = obj.target.counts[2:end] ./ sum(obj.target.counts[2:end])
	sum((overlaps_norm .- target_norm).^2)
end

# on a field of [0,1) x [0,1)
function makebitmapcircle(x, y, r, size)
	bm = falses(size,size)
	pixel = 1/size

	# walk rows from -r to +r, doing the trig to find the number of bits
	# left and right of x to set to true
	for yoffset in -r:pixel:r
		# figure out the range for this row
		alpha = r * sqrt(1 - (yoffset/r)^2) # a big of trig
		#@show alpha
		# convert into bitmap coordinates
		yoffset_bm = iround((y + yoffset) * size + 1)
		#@show yoffset_bm
		if 1 <= yoffset_bm <= size 	# if Y is inside the box
			# convert X into bitmap coordinates, bounding
			xrange_bm = iround(max(1,(x - alpha) * size + 1)) : iround(min(size,(x + alpha) * size + 1))
			#@show xrange_bm
			if (length(xrange_bm) > 0)
				bm[yoffset_bm, xrange_bm] = true
			end
		end
	end
	bm
end

function showbitmap(bm)
	for r in 1:size(bm,1)
		for c in 1:size(bm,2)
			print(bm[r,c] ? '￭' : '·')
		end
		println("")
	end
end

function render(fn, obj::EulerObject, state::EulerState; size=500.0)
	c = CairoSVGSurface(fn, size, size);
	cr = CairoContext(c);
	set_line_width(cr, size / 200.0);
	select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL,
	                     Cairo.FONT_WEIGHT_NORMAL);
	set_font_size(cr, size / 15.0);
	function calcxy(i)
		x = state[(i-1)*2+1]*size
		y = state[(i-1)*2+2]*size
		x,y
	end
	for i = 1:3
		set_source_rgba(cr, i==1?1:0, i==2?1:0, i==3?1:0, .5)
		(x, y) = calcxy(i)
		arc(cr, x, y, obj.sizes[i]*size, 0, 2*pi)
		fill(cr)
	end
	for i = 1:3 # put edges and labels on top
		set_source_rgba(cr, i==1?1:0, i==2?1:0, i==3?1:0, 1)
		(x, y) = calcxy(i)
		arc(cr, x, y, obj.sizes[i]*size, 0, 2*pi)
		stroke(cr);
		extents = text_extents(cr, obj.labels[i]);	
		move_to(cr, x-(extents[3]/2 + extents[1]), y-(extents[4]/2 + extents[2]));
		show_text(cr, obj.labels[i]);
	end
	finish(c)
end

end # module
