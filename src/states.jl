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

# given count_totals, which is unscaled sizes, compute a size which 
# depends on the shape type, along with bounding box info.
# lb/ub for x/y params are half the maximum extent
# lb/ub for other parameters are usually [0,1]
function compute_shape_sizes(specs, count_totals, sizesum)
	@assert length(specs) == length(count_totals)
	shape_sizes = similar(count_totals, Float64)
	bounding_boxes = Float64[] #similar(count_totals, Float64)

	# rescale count_totals (areas) to sum to sizesum
	areas = sizesum * count_totals / sum(count_totals)

	for i in 1:length(specs)
		if specs[i].shape == :circle
			# area = pi r^2; r = sqrt(area/pi)
			shape_sizes[i] = sqrt(areas[i] / pi)
			append!(bounding_boxes, [shape_sizes[i], shape_sizes[i]])
		elseif specs[i].shape == :square
			# area = side^2; halfside = sqrt(area) / 2
			shape_sizes[i] = sqrt(areas[i]) / 2
			append!(bounding_boxes, [shape_sizes[i], shape_sizes[i]])
		elseif specs[i].shape == :triangle
			# area = side^2 * sqrt(3)/4; dist to vertex (bounding circle radius) = side/2 / cos(pi/6)
			side = sqrt(4 * areas[i] / sqrt(3))
			shape_sizes[i] = side / (2cos(pi/6))
			append!(bounding_boxes, [shape_sizes[i], shape_sizes[i]]) # I think side/2 is right, but get errors?
		elseif specs[i].shape == :rectangle
			# for rectangles, one of the parameters is the eccentricity, which
			# varies from 0 to 1. This maps to a height/width ratio that varies
			# from 1:4 to 4:1 (capped so that boundaries can be drawn). The
			# conversion formula is hw = 16 ^ (ecc - 0.5). 
			# From there, height = sqrt(area) * hw
			# and width = sqrt(area) / hw. This means that shape_size for a 
			# rectangle is sqrt(area)/2.
			shape_sizes[i] = sqrt(areas[i])/2
			append!(bounding_boxes, [2shape_sizes[i], 2shape_sizes[i], 0])
				# to give room depending on ecc; no need to bound eccentricity
		else
			error("unknown shape: ", specs[i].shape)
		end
	end
	shape_sizes, bounding_boxes
end

"""
    make_euler_object(labels, counts, specs; sizesum=1) -> EulerObject

Configure an `EulerObject` with the data to plot and the desired style.  For N
classes, `labels` is a vector of length N specifying their names; `counts` is a
boolean matrix with N columns for which each row is an observation; and `specs`
is an `EulerSpec` or N-vector thereof with the shape.  The keyword argument
`sizenum` controls the areas of the shapes.
"""
function make_euler_object(labels, counts, specs::Vector{EulerSpec}; sizesum = 1)
	target = DisjointSet(counts, labels)
	count_totals = vec(sum(counts, dims=1))

	@assert length(labels) == length(count_totals)

	# use the spec to figure out how to translate to and from a state
	nparams = sum([length(spec.clamp) for spec in specs])
	
	update_statepos!(specs) # this lets us look the up parameters in the state vector
	
	# we have non-normalized areas, but need to compute a per-shape size (e.g., radius) 
	# that indicates the maximum distance from the object center to the farthest edge
	shape_sizes, bounding_boxes = compute_shape_sizes(specs, count_totals, sizesum)
	# then, use those sizes to generate bounds

	# this forces the centers to not overlap with the boundary
	# and inserts any non-NaN clamps into the vector
	clamp_vec = vcat([sp.clamp for sp in specs]...)
	#@show clamp_vec
	#@show shape_sizes
	lb = ifelse.(isnan.(clamp_vec), bounding_boxes, clamp_vec)
	ub = ifelse.(isnan.(clamp_vec), 1 .- bounding_boxes, clamp_vec)
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
	make_euler_object(labels, counts, EulerSpec[deepcopy(spec)::EulerSpec for x in 1:length(labels)]; q...)

"""
    random_state(eo::EulerObject)

Pick random parameters within the upper and lower bounds of each `EulerSpec` in `eo`.
"""
random_state(eo::EulerObject) = rand(eo.nparams) .* (eo.ub .- eo.lb) .+ eo.lb
