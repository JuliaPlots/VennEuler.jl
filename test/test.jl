using VennEuler, Base.Test

# just use @test macro
# run with:
# reload("src/VennEuler.jl"); include("test/test.jl")

srand(2)

# set up an initial test of 3 circles with random set membership

# data is a boolean matrix of N rows/observations, with M columns/sets,
# and any number of set memberships for each observation (although often one
# or more in practice). 
randdata = randbool(20, 3) # 3 cols
setlabels = ["A", "B", "C"]
# specification is either a single object (if all sets are treated the same)
# or a list of objects (to treat them differently)
spec = EulerSpec(:circle, [NaN, NaN])
specdefault = EulerSpec()
@test isequal(spec.shape, specdefault.shape)
@test isequal(spec.clamp, specdefault.clamp)
@test_throws ErrorException EulerSpec(:tesseract)
@test_throws ErrorException EulerSpec(:circle, [1, 1, 1])

# make sure DisjointSets can be constructed
ds1 = DisjointSet(randdata, setlabels)
@test ds1.counts == [0,0,5,3,5,1,4,2]

# make an Euler object
# es, eo = makeeulerobject(setlabels, vec(sum(randdata,1)), ds1, sizesum=.5) # old way
# @test_approx_eq_eps(eo.lb, [0.17939, 0.17939, 0.193763, 0.193763, 0.126848, 0.126848], .001)
# @test all(eo.lb .<= es .<= eo.ub)
es, eo = makeeulerobject(setlabels, randdata, sizesum=.5) # test shortcut
@test_approx_eq_eps(eo.lb, [0.17939, 0.17939, 0.193763, 0.193763, 0.126848, 0.126848], .001)
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

es2, eo2 = makeeulerobject(simplelabels, simpledata, sizesum=.5)
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


