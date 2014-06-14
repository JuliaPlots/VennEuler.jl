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
