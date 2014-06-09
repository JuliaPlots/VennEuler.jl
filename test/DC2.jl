using VennEuler

data, labels = readcsv("test/DC2.csv", header=true)
data = bool(data)
labels = vec(labels)

areas = vec(sum(data,1))
radii = sqrt(areas / pi) # a = pi r^2; r = sqrt(a/pi)
es, eo = makeeulerobject(labels, radii, DisjointSet(data, labels))

(minf,minx,ret) = optimize(eo, es, ftol=-1, xtol=0.0025, maxtime=240, pop=2000)
println("got $minf at $minx (returned $ret)")

render("test/DC2.svg", eo, minx)
