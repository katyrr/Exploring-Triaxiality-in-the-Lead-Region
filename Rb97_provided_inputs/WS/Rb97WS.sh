# Test bash script to automate WS programs for Rb 97
# Run from directory /code/Rb97/WS/Outputs


echo running WS programs for Rb 97

./../../../Executables/WS/swgamma < ../Inputs/Rb97SWG.DAT

echo swgamma done

./../../../Executables/WS/wsdcup <  ../Inputs/Rb97WSD.DAT

echo wsdcup done

./../../../Executables/WS/asyrws < ../Inputs/Rb97ASY.DAT

echo asyrws done

./../../../Executables/WS/probaws < ../Inputs/Rb97PR.DAT

echo probaws done