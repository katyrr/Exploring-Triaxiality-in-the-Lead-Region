echo running ASYRMO ...
./../../../Executables/MO/asyrmo < ../Inputs/ASY_Ir191_eps_-0.300_gamma_30.DAT
cp ASYRMO.out ASY_Ir191_eps_-0.300_gamma_30.OUT
./../../../Executables/MO/asyrmo < ../Inputs/ASY_Ir191_eps_-0.100_gamma_30.DAT
cp ASYRMO.out ASY_Ir191_eps_-0.100_gamma_30.OUT
./../../../Executables/MO/asyrmo < ../Inputs/ASY_Ir191_eps_0.100_gamma_30.DAT
cp ASYRMO.out ASY_Ir191_eps_0.100_gamma_30.OUT
./../../../Executables/MO/asyrmo < ../Inputs/ASY_Ir191_eps_0.300_gamma_30.DAT
cp ASYRMO.out ASY_Ir191_eps_0.300_gamma_30.OUT

echo done