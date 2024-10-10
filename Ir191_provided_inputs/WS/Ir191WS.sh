# Test bash script to automate WS programs for Ir 191
# Run from directory /code/Ir191/WS/Outputs


echo running WS programs for Ir 191

./../../../Executables/WS/swgamma < ../Inputs/IR191SWG.DAT

echo swgamma done

./../../../Executables/WS/wsdcup <  ../Inputs/IR191WSD.DAT

echo wsdcup done

./../../../Executables/WS/asyrws < ../Inputs/IR191ASY.DAT

echo asyrws done

./../../../Executables/WS/probaws < ../Inputs/IR191PR.DAT

echo probaws done