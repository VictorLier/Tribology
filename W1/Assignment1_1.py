
import math
from math import sqrt

# film parameters
Lambda_hydro=5
Lambda_partial=1

# Ra parameters
Ra_mill_min = 0.8
Ra_mill_max = 6.3

Ra_grind_min = 0.1
Ra_grind_max = 0.4

# Rq parameters
Rq_mill_min = Ra_mill_min * 1.25
Rq_mill_max = Ra_mill_max * 1.25

Rq_grind_min = Ra_grind_min * 1.25
Rq_grind_max = Ra_grind_max * 1.25

# Film thickness parameters
h_min_mill_hydro = sqrt(2) * Lambda_hydro * Rq_mill_min
h_max_mill_hydro = sqrt(2) * Lambda_hydro * Rq_mill_max

h_min_mill_partial = sqrt(2) * Lambda_partial * Rq_mill_min
h_max_mill_partial = sqrt(2) * Lambda_partial * Rq_mill_max

h_min_grind_hydro = sqrt(2) * Lambda_hydro * Rq_grind_min
h_max_grind_hydro = sqrt(2) * Lambda_hydro * Rq_grind_max

h_min_grind_partial = sqrt(2) * Lambda_partial * Rq_grind_min
h_max_grind_partial = sqrt(2) * Lambda_partial * Rq_grind_max

# print results
print(f"Film thickness for mill finish with hydrodynamic lubrication: {h_min_mill_hydro:.3f} - {h_max_mill_hydro:.3f}")
print(f"Film thickness for mill finish with partial lubrication: {h_min_mill_partial:.3f} - {h_max_mill_partial:.3f}")
print(f"Film thickness for ground finish with hydrodynamic lubrication: {h_min_grind_hydro:.3f} - {h_max_grind_hydro:.3f}")
print(f"Film thickness for ground finish with partial lubrication: {h_min_grind_partial:.3f} - {h_max_grind_partial:.3f}")

