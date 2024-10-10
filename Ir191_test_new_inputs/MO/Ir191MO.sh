# Test bash script to automate MO programs for Ir191
# Run from directory /code/Ir191/MO/Outputs


echo running MO programs for Ir 191

./../../../Executables/MO/gampn <  ../Inputs/GAM191IR.DAT

echo gampn done

./../../../Executables/MO/asyrmo < ../Inputs/ASY191IR.DAT

echo asyrmo done

./../../../Executables/MO/probamo < ../Inputs/PROB191IR.DAT

echo probamo done