using VennEuler

data, labels = readcsv("test/DC2.csv", header=true)
data = bool(data)
labels = vec(labels)

eo = make_euler_object(labels, data, EulerSpec()) # circles, for now

(minf,minx,ret) = optimize(eo, random_state(eo), ftol=-1, xtol=0.0025, maxtime=120, pop=1000)
println("got $minf at $minx (returned $ret)")

render("test/DC2.svg", eo, minx)

eo = make_euler_object(labels, data, [EulerSpec(:circle), EulerSpec(:square), EulerSpec(:triangle),
	EulerSpec(:circle), EulerSpec(:square), EulerSpec(:triangle)]) # mixed shapes

(minf,minx,ret) = optimize(eo, random_state(eo), ftol=-1, xtol=0.0025, maxtime=240, pop=1000)
println("got $minf at $minx (returned $ret)")

render("test/DC2-mixed.svg", eo, minx)
