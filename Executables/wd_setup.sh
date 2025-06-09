echo checking directory setup...
if [ ! -x "./64bit/MO/GAMPN.exe" ]; then
chmod +x ./64bit/MO/GAMPN.exe
echo modified permissions to make GAMPN.exe executable
fi
if [ ! -x "./64bit/MO/ASYRMO.exe" ]; then
chmod +x ./64bit/MO/ASYRMO.exe
echo modified permissions to make ASYRMO.exe executable
fi
if [ ! -x "./64bit/MO/PROBAMO.exe" ]; then
chmod +x ./64bit/MO/PROBAMO.exe
echo modified permissions to make PROBAMO.exe executable
fi
if [ ! -d "../Pt177/Inputs" ]; then
mkdir ../Pt177/Inputs
fi
if [ ! -d "../Pt177/Scripts" ]; then
mkdir ../Pt177/Scripts
fi
if [ ! -d "../Pt177/Run" ]; then
mkdir ../Pt177/Run
fi
if [ ! -d "../Pt177/Outputs" ]; then
mkdir ../Pt177/Outputs
fi
if [ ! -d "../Pt177/Run/Batch1" ]; then
mkdir ../Pt177/Run/Batch1
fi
if [ ! -d "../Pt177/Run/Batch2" ]; then
mkdir ../Pt177/Run/Batch2
fi
if [ ! -d "../Pt177/Run/Batch3" ]; then
mkdir ../Pt177/Run/Batch3
fi
if [ ! -d "../Pt177/Run/Batch4" ]; then
mkdir ../Pt177/Run/Batch4
fi
if [ ! -d "../Pt177/Run/Batch5" ]; then
mkdir ../Pt177/Run/Batch5
fi
if [ ! -d "../Pt177/Run/Batch6" ]; then
mkdir ../Pt177/Run/Batch6
fi
if [ ! -d "../Pt177/Run/Batch7" ]; then
mkdir ../Pt177/Run/Batch7
fi
if [ ! -d "../Pt177/Run/Batch8" ]; then
mkdir ../Pt177/Run/Batch8
fi