using VennEuler

data, labels = readcsv("test/DC2.csv", header=true)
data = bool(data)
labels = vec(labels)

eo = make_euler_object(labels, data, EulerSpec()) # circles, for now

es = random_state(eo)
VennEuler.eval_euler_state(eo, es)
@profile for i in 1:1000 VennEuler.eval_euler_state(eo, es) end

Profile.print()