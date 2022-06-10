x = collect(0:10)
y = map(exp, 0:10)
m, c = SomaticEvolution.fitexponential(x, y)
@test m ≈ 1
@test c == 0

y = map(x -> exp(x), 2:2:22)
m, c = SomaticEvolution.fitexponential(x, y, false)
@test m ≈ 2
@test c ≈ 2