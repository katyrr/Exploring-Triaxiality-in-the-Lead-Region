echo running GAMPN ...
./../../../Executables/MO/gampn < ../Inputs/GAM_Ir191_eps_-0.300_gamma_30.DAT
cp GAMPN.out GAM_Ir191_eps_-0.300_gamma_30.OUT
./../../../Executables/MO/gampn < ../Inputs/GAM_Ir191_eps_-0.100_gamma_30.DAT
cp GAMPN.out GAM_Ir191_eps_-0.100_gamma_30.OUT
./../../../Executables/MO/gampn < ../Inputs/GAM_Ir191_eps_0.100_gamma_30.DAT
cp GAMPN.out GAM_Ir191_eps_0.100_gamma_30.OUT
./../../../Executables/MO/gampn < ../Inputs/GAM_Ir191_eps_0.300_gamma_30.DAT
cp GAMPN.out GAM_Ir191_eps_0.300_gamma_30.OUT

echo done