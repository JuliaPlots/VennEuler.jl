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
	render,
	random_state


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
		if (in(shape, [:circle, :square, :triangle]))
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
	if in(shape, [:circle, :square, :triangle])
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

function compute_shape_sizes(specs, count_totals, sizesum)
	@assert length(specs) == length(count_totals)
	shape_sizes = similar(count_totals, Float64)
	for i in 1:length(specs)
		if specs[i].shape == :circle
			# area = pi r^2; r = sqrt(area/pi)
			shape_sizes[i] = sqrt(count_totals[i] / pi) 
		elseif specs[i].shape == :square
			# area = side^2; halfside = sqrt(area) / 2
			shape_sizes[i] = sqrt(count_totals[i]) / 2
		elseif specs[i].shape == :triangle
			# area = side^2 * sqrt(3)/4; dist to vertex = side/2 / cos(pi/6)
			side = sqrt(4 * count_totals[i] / sqrt(3))
			shape_sizes[i] = side / (2 * cos(pi/6))
		else
			error("unknown shape: ", specs[i].shape)
		end
	end
	sizesum * shape_sizes / sum(shape_sizes)
end

function make_euler_object(labels, counts, specs::Vector{EulerSpec}; sizesum = 1)
	target = DisjointSet(counts, labels)
	count_totals = vec(sum(counts,1))

	@assert length(labels) == length(count_totals)

	# use the spec to figure out how to translate to and from a state
	nparams = sum([length(spec.clamp) for spec in specs])
	
	update_statepos!(specs) # this lets us look the up parameters in the state vector
	
	# we have non-normalized areas, but need to compute a per-shape size (e.g., radius) 
	# that indicates the maximum distance from the object center to the farthest edge
	shape_sizes = compute_shape_sizes(specs, count_totals, sizesum)
	# then, use those sizes to generate bounds

	# this forces the centers to not overlap with the boundary
	# and inserts any non-NaN specs into the 
	clamp_vec = [[sp.clamp for sp in specs]...]
	lb = ifelse(isnan(clamp_vec), dupeelements(shape_sizes), clamp_vec)
	ub = ifelse(isnan(clamp_vec), dupeelements(1 .- shape_sizes), clamp_vec)
	@assert all(lb .<= ub)

	# return: state vector, state object (with bounds, closure, etc)
	eo = EulerObject(nparams, labels, lb, ub, shape_sizes, target, specs, identity)
	eo.evalfn = (x,g) -> begin
		@assert length(g) == 0
		cost = VennEuler.eval_euler_state(eo, x) 
		#println(x, " (", cost, ")")
		cost
	end
	eo
end
make_euler_object(labels, counts, spec::EulerSpec; q...) = 
	make_euler_object(labels, counts, [deepcopy(x)::EulerSpec for x in repeated(spec, length(labels))]; q...)

function random_state(eo::EulerObject)
	rand(eo.nparams) .* (eo.ub .- eo.lb) .+ eo.lb
end

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
	# generate a 2-D bitmap from each set
	# foreach element of the DisjointSet, calculate the size of the overlap of the bitmaps
	# compare the overlaps with the target, returning the error metric

	# draws different shapes depending on spec
	bitmaps = [make_bitmap(state[obj.specs[i].statepos], obj.sizes[i], obj.specs[i], px) 
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
function make_bitmap(params, size, spec, px)
	if spec.shape == :circle
		make_bitmap_circle(params[1], params[2], size, px)
	elseif spec.shape == :square
		make_bitmap_square(params[1], params[2], size, px)
	elseif spec.shape == :triangle
		make_bitmap_triangle(params[1], params[2], size, px)
	else
		error("unknown shape: ", spec.shape)
	end
end

function make_bitmap_triangle(x, y, r, px) # r = radius of circle around the triangle
	bm = falses(px,px)
	pixel = 1/px

	# this shape isn't symmetric -- need to be careful that 0,0 is UPPER-LEFT corner, like Cairo

	# TODO: might be able to share some of this with the Cairo rendering

	# walk columns from left to right, doing the trig to find the number of bits
	# to set to true
	h = cos(pi/6) * r # half the side of the triangle
	j = sin(pi/6) * r # distance from center to side 
	for xoffset in -h:pixel:h
		# figure out the range for this column
		p = h - abs(xoffset) # dist from the corner to this column, mirrored
		yoffset = (j - (r+j) * p/h, j) # congruent triangles give ratio of vertical offset from j

		# covert into bitmap coords
		xoffset_bm =  iround((x + xoffset) * px + 1)
		yrange_bm = iround(max(1, (y + yoffset[1]) * px + 1)) : iround(min(px, (y + yoffset[2]) * px + 1))
		bm[yrange_bm, xoffset_bm] = true
	end
	bm
end

function make_bitmap_square(x, y, halfsz, px)
	bm = falses(px, px) # DRY...
	pixel = 1/px

	# compute column and row ranges, then fast assignment
	xrange_bm = iround(max(1, (x - halfsz) * px + 1)) : iround(min(px, (x + halfsz) * px + 1))
	yrange_bm = iround(max(1, (y - halfsz) * px + 1)) : iround(min(px, (y + halfsz) * px + 1))
	bm[yrange_bm, xrange_bm] = true # matrices are indexed Y,X...
	bm
end

function make_bitmap_circle(x, y, r, px)
	bm = falses(px,px)
	pixel = 1/px

	# walk rows from -r to +r, doing the trig to find the number of bits
	# left and right of x to set to true
	for yoffset in -r:pixel:r
		# figure out the range for this row
		alpha = r * sqrt(1 - (yoffset/r)^2) # a big of trig
		#@show alpha
		# convert into bitmap coordinates
		yoffset_bm = iround((y + yoffset) * px + 1)
		#@show yoffset_bm
		if 1 <= yoffset_bm <= px 	# if Y is inside the box
			# convert X into bitmap coordinates, bounding
			xrange_bm = iround(max(1,(x - alpha) * px + 1)) : iround(min(px,(x + alpha) * px + 1))
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

get_params(state::EulerState, spec::EulerSpec) = tuple(state[spec.statepos]...)

reds =   [1,  0, 0, .7, 0, .7, .2, .7]
greens = [0, .8, 0, 0, .7, .7, .2, .7]
blues =  [0,  0, 1, .7, .7, 0, .2, .7]	

function setup_shape!(cr, spec, state, size, px) 
	if spec.shape == :circle
		x,y = get_params(state, spec)
		move_to(cr, x*px+size*px, y*px)
		arc(cr, x*px, y*px, size*px, 0, 2*pi)
	elseif spec.shape == :square
		x,y = get_params(state, spec)
		halfside = size
		rectangle(cr, (x - halfside)*px, (y - halfside)*px, 2halfside*px, 2halfside*px)
	elseif spec.shape == :triangle
		x,y = get_params(state, spec)
		# convert into our triangle notation
		r = size
		h = cos(pi/6) * r # half the side of the triangle
		j = sin(pi/6) * r # distance from center to side 
		# start at lower-left, then top, then lower-right, then close
		move_to(cr, (x-h)*px, (y+j)*px)
		line_to(cr, x*px, (y-r)*px)
		line_to(cr, (x+h)*px, (y+j)*px)
		close_path(cr)
	end
	x,y
end

function render(fn, obj::EulerObject, state::EulerState; verbose=0, px=500.0)
	c = CairoSVGSurface(fn, px, px);
	cr = CairoContext(c);
	set_line_width(cr, px / 200.0);
	select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL,
	                     Cairo.FONT_WEIGHT_NORMAL);
	set_font_size(cr, px / 25.0);
	
	@assert length(obj.labels) <= length(reds)

	# foreach object, set the color, calc the center, then draw the appropriate object
	for i = 1:length(obj.labels)
		if verbose>0 @show obj.specs[i] end

		set_source_rgba(cr, reds[i], greens[i], blues[i], .3)
		setup_shape!(cr, obj.specs[i], state, obj.sizes[i], px)
		fill(cr)
	end
	for i = 1:length(obj.labels) # put edges and labels on top
		set_source_rgba(cr, reds[i], greens[i], blues[i], 1)

		# could use set_operator(CAIRO_OPERATOR_DEST_OVER) to draw underneath existing stuff
		x,y = setup_shape!(cr, obj.specs[i], state, obj.sizes[i], px)
		stroke(cr);

		extents = text_extents(cr, obj.labels[i]);	
		move_to(cr, x*px-(extents[3]/2 + extents[1]), y*px-(extents[4]/2 + extents[2]));
		show_text(cr, obj.labels[i]);
	end
	finish(c)
end

end # module
