using VennEuler

data, labels = readcsv("test/DC2.csv", header=true)
data = bool(data)
labels = vec(labels)

es, eo = make_euler_object(labels, data, EulerSpec()) # circles, for now

(minf,minx,ret) = optimize(eo, es, ftol=-1, xtol=0.0025, maxtime=240, pop=2000)
println("got $minf at $minx (returned $ret)")

render("test/DC2.svg", eo, minx)
