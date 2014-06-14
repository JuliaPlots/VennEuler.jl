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

