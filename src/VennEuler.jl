module VennEuler

using Base
using Iterators
using NLopt
using Cairo

export
	EulerState,
	EulerObject,
	EulerSpec,
	make_euler_object,
	optimize,
	render


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
		if (shape == :circle)
			@assert length(clamp) == 2
			@assert length(statepos) == 2
		else 
			error("Unknown EulerSpec shape: ", string(shape))
		end
		new(shape, clamp, statepos)
	end
end

EulerSpec() = EulerSpec(:circle)
function EulerSpec(shape)
	if (shape == :circle)
		EulerSpec(shape, [NaN, NaN], [0, 0])
	else
		error("Unknown EulerSpec shape: ", string(shape))
	end
end

dupeelements(qq) = [qq[ceil(i)] for i in (1:(2*length(qq)))/2]

function update_statepos!(specs::Vector{EulerSpec})
	i = 1
	for spec in specs
		for j = 1:length(spec.statepos)
			spec.statepos[j] = i
			i += 1
		end
		#@show spec
	end
end

function make_euler_object(labels, counts, specs::Vector{EulerSpec}; sizesum = 1)
	target = DisjointSet(counts, labels)
	sizes = vec(sum(counts,1))

	@assert length(labels) == length(sizes)

	# use the spec to figure out how to translate to and from a state
	nparams = sum([length(spec.clamp) for spec in specs])
	# TODO: turn the below into 
	
	
	# convert areas to radii, then normalize
	radii = sqrt(sizes / pi)
	radii = sizesum * radii / sum(radii)
	# this forces the centers to not overlap with the boundary
	lb = dupeelements(radii)
	ub = dupeelements(1 .- radii)
	@assert all(lb .< ub)
	# create the initial state vector
	es::EulerState = rand(nparams) .* (ub .- lb) .+ lb
	
	# return: state vector, state object (with bounds, closure, etc)
	eo = EulerObject(nparams, labels, lb, ub, radii, target, specs, identity)
	eo.evalfn = (x,g) -> begin
		@assert length(g) == 0
		cost = VennEuler.eval_euler_state(eo, x) 
		#println(x, " (", cost, ")")
		cost
	end
	(es, eo)
end
make_euler_object(labels, counts, spec::EulerSpec; q...) = 
	make_euler_object(labels, counts, [deepcopy(x)::EulerSpec for x in repeated(spec, length(labels))]; q...)

function optimize(obj::EulerObject, state::EulerState; xtol=1/200, ftol=1.0e-7, maxtime=30, init_step=.1,
	alg = :GN_CRS2_LM, pop = 0)
	opt = Opt(alg, length(state))
	lower_bounds!(opt, obj.lb)
	upper_bounds!(opt, obj.ub)
	xtol_abs!(opt, xtol)
	ftol_abs!(opt, ftol)
	min_objective!(opt, obj.evalfn)
	maxtime!(opt, maxtime)
	initial_step!(opt, init_step)
	population!(opt, pop)
	NLopt.optimize(opt, state)
end

function eval_euler_state(obj::EulerObject, state::EulerState; verbose::Int64=0, px::Int64=200)
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
	set_font_size(cr, size / 25.0);
	function calcxy(i)
		x = state[(i-1)*2+1]*size
		y = state[(i-1)*2+2]*size
		x,y
	end
	reds =   [1,  0, 0, .7, 0, .7, .2, .7]
	greens = [0, .8, 0, 0, .7, .7, .2, .7]
	blues =  [0,  0, 1, .7, .7, 0, .2, .7]
	@assert length(obj.labels) <= length(reds)
	for i = 1:length(obj.labels)
		set_source_rgba(cr, reds[i], greens[i], blues[i], .3)
		(x, y) = calcxy(i)
		arc(cr, x, y, obj.sizes[i]*size, 0, 2*pi)
		fill(cr)
	end
	for i = 1:length(obj.labels) # put edges and labels on top
		set_source_rgba(cr, reds[i], greens[i], blues[i], 1)
		(x, y) = calcxy(i)
		move_to(cr, x+obj.sizes[i]*size, y)
		arc(cr, x, y, obj.sizes[i]*size, 0, 2*pi)
		stroke(cr);
		extents = text_extents(cr, obj.labels[i]);	
		move_to(cr, x-(extents[3]/2 + extents[1]), y-(extents[4]/2 + extents[2]));
		show_text(cr, obj.labels[i]);
	end
	finish(c)
end

end # module
