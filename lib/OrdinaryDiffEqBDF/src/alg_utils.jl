isadaptive(alg::DABDF2) = true
isadaptive(alg::DFBDF) = true

alg_extrapolates(alg::ABDF2) = true
alg_extrapolates(alg::SBDF) = true
alg_extrapolates(alg::MEBDF2) = true
alg_extrapolates(alg::DABDF2) = true

get_current_adaptive_order(alg::FBDF, cache) = cache.order
get_current_alg_order(alg::FBDF, cache) = cache.order

alg_order(alg::TRBDF2) = 2
alg_order(alg::ABDF2) = 2
alg_order(alg::FBDF) = 1 #dummy value
alg_order(alg::SBDF) = alg.order
alg_order(alg::MEBDF2) = 2

alg_order(alg::DABDF2) = 2
alg_order(alg::DFBDF) = 1#dummy value
qsteady_min_default(alg::FBDF) = 9 // 10
qsteady_max_default(alg::FBDF) = 2 // 1

isesdirk(alg::TRBDF2) = true