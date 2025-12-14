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