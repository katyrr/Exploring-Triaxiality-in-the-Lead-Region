echo running PROBAMO ...
./../../../Executables/MO/probamo < ../Inputs/PROB_Ir191_eps_-0.300_gamma_30.DAT
cp PROBAMO.out PROB_Ir191_eps_-0.300_gamma_30.OUT
./../../../Executables/MO/probamo < ../Inputs/PROB_Ir191_eps_-0.100_gamma_30.DAT
cp PROBAMO.out PROB_Ir191_eps_-0.100_gamma_30.OUT
./../../../Executables/MO/probamo < ../Inputs/PROB_Ir191_eps_0.100_gamma_30.DAT
cp PROBAMO.out PROB_Ir191_eps_0.100_gamma_30.OUT
./../../../Executables/MO/probamo < ../Inputs/PROB_Ir191_eps_0.300_gamma_30.DAT
cp PROBAMO.out PROB_Ir191_eps_0.300_gamma_30.OUT

echo done