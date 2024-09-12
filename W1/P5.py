
# In a dry rubbing bearing, a steel shaft is paired with a PTFE sleeve.
# Using the expression by Gohar et.al. p. 65/66 give an estimate of the adhesive coefficient of friction.
# Compare this to the ploughing coefficient of friction if it is assumed that the asperity 
# radius is 100μm and the groove depth is h = 1μm. How does this compare to data sheet values?
# A PTFE reference could be: https://www.bearingworks.com/uploaded-assets/pdfs/retainers/ptfe-datasheet.pdf

# adhesive coefficient of friction:
k = 5e6 # [Pa] - shear strength of PTFE
H = 15e69 # [Pa] - compressive yield strength of PTFE
alpha = (H/k)**2
c = 0.2 # [-] - Assuming some surface contamination (clean: c = 1)
mu_a1 = c/((alpha*(1-c**2))**0.5)
mu_a2 = c*k/H
print(f"Adhesive coefficient of friction: {mu_a1:.2f} or {mu_a2:.2f}")

# ploughing coefficient of friction:
h = 1e-6
R = 100e-6
mu_p = 0.6*(h/R)**0.5

print(f"Ploughing coefficient of friction: {mu_p:.2f}")

print(f"This compares well to the data sheet that states a coefficient of friction of 0.06-0.1")