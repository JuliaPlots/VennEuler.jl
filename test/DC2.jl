using VennEuler

data, labels = readcsv("test/DC2.csv", header=true)
data = bool(data)
labels = vec(labels)

eo = make_euler_object(labels, data, EulerSpec()) # circles, for now

(minf,minx,ret) = optimize(eo, random_state(eo), ftol=-1, xtol=0.0025, maxtime=120, pop=1000)
println("got $minf at $minx (returned $ret)")

render("test/DC2.svg", eo, minx)

# do it again, but with different shapes, and with DSDC (the 2nd one) clamped in the center
eo = make_euler_object(labels, data, [EulerSpec(:circle), EulerSpec(:square, [.5, .5], [0, 0]), 
	EulerSpec(:triangle),
	EulerSpec(:circle), EulerSpec(:square), EulerSpec(:triangle)]) # mixed shapes

(minf,minx,ret) = optimize(eo, random_state(eo), ftol=-1, xtol=0.0025, maxtime=240, pop=1000)
println("got $minf at $minx (returned $ret)")

render("test/DC2-mixed.svg", eo, minx)

# do it again, but with rectangles (DSDC fixed)
eo = make_euler_object(labels, data, [EulerSpec(:rectangle), EulerSpec(:rectangle, [.5, .5, .4], [0, 0, 0]), 
	EulerSpec(:rectangle),
	EulerSpec(:rectangle), EulerSpec(:rectangle), EulerSpec(:rectangle)],
	sizesum=.3)

#(minf,minx,ret) = optimize(eo, random_state(eo), ftol=-1, xtol=0.0025, maxtime=240, pop=1000)
(minf,minx,ret) = optimize_iteratively(eo, random_state(eo), ftol=-1, xtol=0.0025, maxtime=5, pop=100)
println("phase 1: got $minf at $minx (returned $ret)")
(minf,minx,ret) = optimize(eo, minx, ftol=-1, xtol=0.001, maxtime=30, pop=100)
println("phase 2: got $minf at $minx (returned $ret)")

render("test/DC2-rects.svg", eo, minx)
