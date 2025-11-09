# Session Notes – CEF vs Effective Doublet (2025-02-15)

- Verified that the current `MF_Er_CaWO4_v1b.m` CEF branch (lines 105-205) is mathematically self-consistent: `cf.m` builds the Stevens Hamiltonian and `Ising_proj.m` diagonalises `H_cf + g_J \mu_B \mathbf{B}\cdot \mathbf{J}` before projecting the lowest Kramers pair.
- Projecting `J_{x,y,z}` onto the zero-field doublet for the literature B-parameters (`[567.02,-945.07,1008.1,0,-22.16,936.29,0.854]`) produces an effective g-tensor of ≈[8.1, 8.1, 4.0]. The effective spin-1/2 model we are fitting against uses [8.3, 8.3, 1.26], so even a perfect fit cannot reconcile the two spectra; the longitudinal component is intrinsically too large.
- The hyperfine tensor embedded in `MF_Er_CaWO4_v1b.m` is already the anisotropic effective tensor (−871, −871, 130 MHz). Feeding it into the full-J Hamiltonian double-counts the anisotropy. The CEF route should start from the microscopic scalar hyperfine constant and let the projection generate the anisotropy.
- Quadrupole interaction in the CEF branch is currently coded as `Q_x I_x^2 + Q_y I_y^2 + Q_z I_z^2`. A realistic electric quadrupole term should be expressed in Stevens form (e.g., `B_2^0 O_2^0 + B_2^2 O_2^2`), otherwise comparisons with measured lineshapes are inconsistent.
- The added scalar `ion.h4 = 5.26e-3` simply shifts all CEF levels; if the intention was to include an `O_4^0` correction from Neda's thesis it should enter as an operator, otherwise it can be dropped.

_Action items_: refit the microscopic CEF parameters so that the projected doublet reproduces the experimental g-tensor; adjust the hyperfine and quadrupole terms accordingly before attempting further ZEFOZ matching.
