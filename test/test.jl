using VennEuler, Base.Test

# just use @test macro
# run with:
# reload("src/VennEuler.jl"); include("test/test.jl")

srand(1)
randdata = randbool(20, 3) # 3 cols
setlabels = ["A", "B", "C"]

# convert the random data into a DisjointSet
ds1 = DisjointSet(randdata, setlabels)
@test ds1.counts == [4,3,1,3,3,1,1,4]

# make an Euler object
es, eo = makeeulerobject(setlabels, vec(sum(randdata,1)), ds1)
@test_approx_eq_eps(eo.lb, [0.155172,0.155172,0.155172,0.155172,0.189655,0.189655], .0001)
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

es2, eo2 = makeeulerobject(simplelabels, vec(sum(simpledata,1)), DisjointSet(simpledata, simplelabels))
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


