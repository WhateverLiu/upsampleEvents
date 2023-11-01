

Rcpp::sourceCpp("tests/GA.cpp", verbose = 1)
initxy = runif(2, -10, 10)
tmp = testGA(initxy = initxy, initNoise = 1, minNoise = 1e-6, 
             popuSize = 100, survivalSize = 30, maxGen = 1000, 
             Ngen2minNoise = 70, reproduceSelection = "", 
             randomSeed = 42, maxCore = 7000)


