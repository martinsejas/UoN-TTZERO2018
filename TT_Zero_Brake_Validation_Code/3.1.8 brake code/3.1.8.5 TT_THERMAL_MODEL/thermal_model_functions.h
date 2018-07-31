#ifndef THERMAL_MODEL_FUNCTIONS_H_INCLUDED
#define THERMAL_MODEL_FUNCTIONS_H_INCLUDED

using namespace std;


//Function to read text file and put values on a double vector
void vfReadTextIntoVector(string* strFilePath, vector <double> *vDataVector, int* iIterator);

//Function to print final values into excel
void vfPrintValuesToExcel(int* iCorners, vector <double> *vCornerVelocityLosses, vector <double> *CornerDistances, vector <double> *vEnergyLosses,
                          vector <double> *vdHeatGeneratedByCorner, vector <double> *vdTimeonCorner, vector <double> *vCornerInitialVelocity,
vector <double> *vCornerExitVelocity, vector <double> *vdCornerTempRise);


//Function to calculate angular acceleration
double dbfCalcAngularAcceleration (double *dbGearRatio, int iCurrentIndex, vector <double> *vdRPM);

//Function to calculate sliding force from torque
double dbfCalcForce (double *dbGearRatio, int iCurrentIndex, vector <double> *vdTorque, double *dbFrontWheelRadius);

//Function to calculate the temperature rise on corner
double dbfStepHeatConductedOut (double* dbBrakeDiscK, double* dbBrakeDiscThickness, double* dbPadArea, double* dbHeatGenerated);

//Function to calculate force
double dbfCalculateForce (int z, vector <double> *vdVelocity, double* bike_mass, vector <double> *vdDistance, double* ke);

//function to process test data
void vfProcessBruntTestData(vector <double> *Time, vector <double> *BrakeTemp, vector <double> *Speed, string *strFilePath);

//function to print calculated test data
void vfPrintTestData(int *iSize, vector <double> *vTime, vector <double> *vTemp, vector <double> *vSpeed, vector <double> *vBrakeTemp);

void vfPrintFinalTempData(int iSize, vector <double> *vDistance, vector <double> *vTemp, vector <int> *vVelocity);

#endif // THERMAL_MODEL_FUNCTIONS_H_INCLUDED
