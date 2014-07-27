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
spec = EulerSpec(:circle, [NaN, NaN], [0, 0])
specdefault = EulerSpec()
@test isequal(spec.shape, specdefault.shape)
@test isequal(spec.clamp, specdefault.clamp)
@test_throws ErrorException EulerSpec(:tesseract)
@test_throws ErrorException EulerSpec(:circle, [1, 1, 1], [0, 0, 0])

spec2 = EulerSpec(:circle, [.5, NaN], [0, 0])
specs1 = [deepcopy(spec), deepcopy(spec), spec2]
VennEuler.update_statepos!(specs1)
@test isequal(specs1[1].statepos, [1, 2])
@test isequal(specs1[2].statepos, [3, 4])

ss, bb = VennEuler.compute_shape_sizes(specs1, vec(sum(randdata,1)), .5)
@test_approx_eq_eps(ss, [0.244301,0.263876,0.172747], .001)

# make sure DisjointSets can be constructed
ds1 = VennEuler.DisjointSet(randdata, setlabels)
@test ds1.counts == [0,0,5,3,5,1,4,2]

# make an Euler object
eo = make_euler_object(setlabels, randdata, spec, sizesum=.5) # test shortcut
@test_approx_eq_eps(eo.lb, [0.244301, 0.244301, 0.263876, 0.263876, 0.172747, 0.172747], .001)
es = random_state(eo)
@test all(eo.lb .<= es .<= eo.ub)

@test VennEuler.get_params(es, specs1[1]) == tuple(es[1:2]...)

# tests for bitmap operations
bmc1 = VennEuler.make_bitmap_circle(.5, .5, .2, 20)
bmc2 = VennEuler.make_bitmap_circle(.5, .5, .001, 20)
@test sum(bmc1) == 57
println("")
VennEuler.showbitmap(bmc1)
@test sum(bmc2) == 1

spec_sq = EulerSpec(:square, [NaN, NaN], [0, 0])
bmsq1 = VennEuler.make_bitmap([.5, .5], .2, spec_sq, 20)
@test sum(bmsq1) == 81
println("")
VennEuler.showbitmap(bmsq1)

spec_tr = EulerSpec(:triangle, [NaN, NaN], [0, 0])
bmtr1 = VennEuler.make_bitmap([.5, .5], .15, spec_tr, 20)
#@show bmtr1 #sum(bmsq1) == 81
println("")
VennEuler.showbitmap(bmtr1)

spec_rect = EulerSpec(:rectangle, [NaN, NaN, NaN], [0, 0, 0])
bmrect1 = VennEuler.make_bitmap([.5, .5, 1], .2, spec_rect, 20)
bmrect2 = VennEuler.make_bitmap([.5, .5, 0], .2, spec_rect, 20)
println("")
VennEuler.showbitmap(bmrect1)
println("")
VennEuler.showbitmap(bmrect2)

println("")

# evaluating Euler states
simplelabels = ["A", "B"]
simpledata = bool([1 0; 1 1; 0 1])

eo2 = make_euler_object(simplelabels, simpledata, spec, sizesum=.5)
@show VennEuler.eval_euler_state(eo2, [.33, .5, .66, .5])
@show eo2.evalfn([.33, .5, .66, .5], [])
@time VennEuler.eval_euler_state(eo2, [.38, .5, .62, .5])

# optimization
(minf,minx,ret) = optimize(eo2, random_state(eo2))
println("got $minf at $minx (returned $ret)")
#@test ret == :FTOL_REACHED

eo = make_euler_object(setlabels, randdata, [spec, spec_tr, spec_sq], sizesum=.5)
(minf,minx,ret) = optimize(eo, random_state(eo), ftol=1/10000)
println("got $minf at $minx (returned $ret)")
#@test ret == :FTOL_REACHED

spec_fixed = EulerSpec(:circle, [0.5, 0.5], [0, 0])
eo = make_euler_object(setlabels, randdata, [spec_fixed, spec_tr, spec_rect], sizesum=.2)
# TODO: better feedback if sizesum is too big
(minf,minx,ret) = optimize(eo, random_state(eo), ftol=1/10000)
println("got $minf at $minx (returned $ret)")
@test isequal(minx[1:2], spec_fixed.clamp)

(minf,minx,ret) = optimize_iteratively(eo, random_state(eo), ftol=1/10000, verbose=true)
println("got $minf at $minx (returned $ret)")

# output
render("test1.svg", eo, minx, verbose=1)


