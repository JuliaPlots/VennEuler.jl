using VennEuler, Base.Test

# just use @test macro

srand(1)
randdata = randbool(20, 3) # 3 cols
setlabels = ["A", "B", "C"]

# convert the random data into a DisjointSet
ds1 = DisjointSet(randdata, setlabels)
@show ds1

# make an Euler object
eo = makeeulerobject(setlabels, vec(sum(randdata,1)), ds1)
@show eo

# tests for bitmap operations
bmc1 = VennEuler.makebitmapcircle(.5, .5, .2, 20)
bmc2 = VennEuler.makebitmapcircle(.5, .5, .001, 20)
@test sum(bmc1) == 57
VennEuler.showbitmap(bmc1)
@test sum(bmc2) == 1

# evaluating Euler states

# optimization


