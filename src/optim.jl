function optimize_iteratively(obj::EulerObject, state::EulerState; xtol=1/200, ftol=1.0e-7, maxtime=30, init_step=.1,
	alg = :GN_ORIG_DIRECT, pop = 0, verbose = false)
	
	order = sortperm(obj.sizes, rev=true)

	# foreach object, in order
	orig_lb = copy(obj.lb)
	orig_ub = copy(obj.ub)
	minf = 0.0; minx = []; ret = 0.0;
	for obj_idx in order
		if (verbose) println("optimizing #$(obj_idx)") end
		# lock every other objects' bounds to the current value
		obj′ = deepcopy(obj)
		obj′.lb = copy(state); obj′.ub = copy(state)
		thisrange = obj.specs[obj_idx].statepos
		if (verbose) println("preserving range: $thisrange") end
		obj′.lb[thisrange] = orig_lb[thisrange]
		obj′.ub[thisrange] = orig_ub[thisrange]
		if (verbose) println("lb = $(obj′.lb)\nub = $(obj′.ub)") end
		# optimize
		(minf,minx,ret) = optimize(obj′, state, xtol=xtol, ftol=ftol, maxtime=maxtime, 
			init_step=init_step, alg=alg, pop=pop)
		println("got $minf at $minx (returned $ret)")
		state = minx
	end
	minf, state, ret
end

function optimize(obj::EulerObject, state::EulerState; xtol=1/200, ftol=1.0e-7, maxtime=30, init_step=.1,
	alg = :GN_ORIG_DIRECT, pop = 0)
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
	not_bitmaps = map(bm->.!bm, bitmaps)

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
		for bm in 1:length(compare_str)
			# TODO: create and!(a,b) that's like bitarray:956, but in-place, for speed
			overlap = overlap .& (compare_str[bm] == '1' ? bitmaps[bm] : not_bitmaps[bm])
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
	elseif spec.shape == :rectangle
		make_bitmap_rectangle(params[1], params[2], params[3], size, px)
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
		xoffset_bm =  round(Integer, (x + xoffset) * px + 1)
		yrange_bm = round(Integer, max(1, (y + yoffset[1]) * px + 1)) : round(Integer, min(px, (y + yoffset[2]) * px + 1))
		bm[yrange_bm, xoffset_bm] = true
	end
	bm
end

function make_bitmap_rectangle(x, y, ecc, halfsz, px)
	bm = falses(px, px) # DRY...
	pixel = 1/px

	hw = 4 ^ (ecc - 0.5) # height/width scaling (sqrt of ratio)
	# compute column and row ranges, then fast assignment
	xrange_bm = round(Integer, max(1, (x - halfsz/hw) * px + 1)) : round(Integer, min(px, (x + halfsz/hw) * px + 1))
	yrange_bm = round(Integer, max(1, (y - halfsz*hw) * px + 1)) : round(Integer, min(px, (y + halfsz*hw) * px + 1))
	bm[yrange_bm, xrange_bm] = true # matrices are indexed Y,X...
	bm
end

function make_bitmap_square(x, y, halfsz, px)
	make_bitmap_rectangle(x, y, .5, halfsz, px)
	# bm = falses(px, px) # DRY...
	# pixel = 1/px

	# # compute column and row ranges, then fast assignment
	# xrange_bm = round(Integer, max(1, (x - halfsz) * px + 1)) : round(Integer, min(px, (x + halfsz) * px + 1))
	# yrange_bm = round(Integer, max(1, (y - halfsz) * px + 1)) : round(Integer, min(px, (y + halfsz) * px + 1))
	# bm[yrange_bm, xrange_bm] = true # matrices are indexed Y,X...
	# bm
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
		yoffset_bm = round(Integer, (y + yoffset) * px + 1)
		#@show yoffset_bm
		if 1 <= yoffset_bm <= px 	# if Y is inside the box
			# convert X into bitmap coordinates, bounding
			xrange_bm = round(Integer, max(1,(x - alpha) * px + 1)) : round(Integer, min(px,(x + alpha) * px + 1))
			#@show xrange_bm
			if (length(xrange_bm) > 0)
				@inbounds bm[yoffset_bm, xrange_bm] = true
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

