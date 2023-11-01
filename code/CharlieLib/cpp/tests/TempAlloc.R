

Rcpp::sourceCpp("tests/TempAlloc.cpp", verbose = 1)


truth = test(100, 5000)
rst = testBenchmark(100, 5000)






