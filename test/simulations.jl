subclones = [
    SomaticEvolution.CloneTracker(
        1, 1, 1.4437821887595856, [2], 2, 1, 1.0, 9
    ),
    SomaticEvolution.CloneTracker(
        1, 1, 1.7344043713478965, [1, 4], 3, 2, 1.6666666666666667, 90
    )
]
N = 100
@test SomaticEvolution.getclonesize(N, subclones) == [1, 9, 90]

#check birth and death rates are calculated correctly
b, d, selection = 1, 0, [1, 2]
brates, drates = SomaticEvolution.set_branching_birthdeath_rates(b, d, selection) 
@test brates == [1, 2, 3]
@test drates == [0, 0, 0]