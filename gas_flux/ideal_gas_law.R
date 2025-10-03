
R_const <- 0.08205736608096 # L⋅atm⋅K−1⋅mol−1

temp <- 0	# C

pressure <- 1	# atm

volume <- 1 # L

molar_mass <- 44.0095	# g·mol−1

molar_mass <- tibble::tribble(
  ~gas,          ~molar_mass,
  "co2",         44.0095,
  "ch4",         16.04246
)
