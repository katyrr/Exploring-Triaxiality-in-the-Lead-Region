echo checking directory setup...
if [ ! -x "./MacOS/MO/gampn" ]; then
chmod +x ./MacOS/MO/gampn
echo modified permissions to make gampn executable
fi
if [ ! -x "./MacOS/MO/asyrmo" ]; then
chmod +x ./MacOS/MO/asyrmo
echo modified permissions to make asyrmo executable
fi
if [ ! -x "./MacOS/MO/probamo" ]; then
chmod +x ./MacOS/MO/probamo
echo modified permissions to make probamo executable
fi
if [ ! -d "../Examples/Inputs" ]; then
mkdir ../Examples/Inputs
fi
if [ ! -d "../Examples/Scripts" ]; then
mkdir ../Examples/Scripts
fi
if [ ! -d "../Examples/Run" ]; then
mkdir ../Examples/Run
fi
if [ ! -d "../Examples/Outputs" ]; then
mkdir ../Examples/Outputs
fi
if [ ! -d "../Examples/Run/Batch1" ]; then
mkdir ../Examples/Run/Batch1
fi
if [ ! -d "../Examples/Run/Batch2" ]; then
mkdir ../Examples/Run/Batch2
fi
if [ ! -d "../Examples/Run/Batch3" ]; then
mkdir ../Examples/Run/Batch3
fi
if [ ! -d "../Examples/Run/Batch4" ]; then
mkdir ../Examples/Run/Batch4
fi
if [ ! -d "../Examples/Run/Batch5" ]; then
mkdir ../Examples/Run/Batch5
fi
if [ ! -d "../Examples/Run/Batch6" ]; then
mkdir ../Examples/Run/Batch6
fi
if [ ! -d "../Examples/Run/Batch7" ]; then
mkdir ../Examples/Run/Batch7
fi
if [ ! -d "../Examples/Run/Batch8" ]; then
mkdir ../Examples/Run/Batch8
fi