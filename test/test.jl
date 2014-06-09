using VennEuler, Base.Test

# just use @test macro
# run with:
# reload("src/VennEuler.jl"); include("test/test.jl")

srand(2)
randdata = randbool(20, 3) # 3 cols
setlabels = ["A", "B", "C"]

# convert the random data into a DisjointSet
ds1 = DisjointSet(randdata, setlabels)
@test ds1.counts == [0,0,5,3,5,1,4,2]

# make an Euler object
es, eo = makeeulerobject(setlabels, vec(sum(randdata,1)), ds1, sizesum=.5)
@test_approx_eq_eps(eo.lb, [0.1875, 0.1875, 0.21875, 0.21875, 0.09375, 0.09375], .001)
@test all(eo.lb .<= es .<= eo.ub)

# tests for bitmap operations
bmc1 = VennEuler.makebitmapcircle(.5, .5, .2, 20)
bmc2 = VennEuler.makebitmapcircle(.5, .5, .001, 20)
@test sum(bmc1) == 57
VennEuler.showbitmap(bmc1)
@test sum(bmc2) == 1

# evaluating Euler states
simplelabels = ["A", "B"]
simpledata = bool([1 0; 1 1; 0 1])

es2, eo2 = makeeulerobject(simplelabels, vec(sum(simpledata,1)), DisjointSet(simpledata, simplelabels),
	sizesum=.5)
@show VennEuler.evaleulerstate(eo2, [.33, .5, .66, .5])
@show eo2.evalfn([.33, .5, .66, .5], [])
@time VennEuler.evaleulerstate(eo2, [.38, .5, .62, .5])

# optimization
(minf,minx,ret) = optimize(eo2, es2)
println("got $minf at $minx (returned $ret)")
@test ret == :FTOL_REACHED

(minf,minx,ret) = optimize(eo, es, ftol=1/10000)
println("got $minf at $minx (returned $ret)")
@test ret == :FTOL_REACHED

# output
render("test1.svg", eo, minx)


