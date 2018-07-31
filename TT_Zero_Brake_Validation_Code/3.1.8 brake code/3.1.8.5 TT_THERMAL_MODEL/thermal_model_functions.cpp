#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include "thermal_model_functions.h"

using namespace std;

void vfReadTextIntoVector(string* strFilePath, vector <double> *vDataVector, int* iIterator)
{

    fstream fData; // declare fstream
    vector <string> vParsingVector; // make a string vector for parsing
    fData.open(strFilePath->c_str()); //open file
    if(fData.good()) // Check if file is open
    {
        while(!fData.eof()) // Loop while we have not reached the end of the file
        {
            string s; //Declare string to read file
            getline(fData, s, ','); // Read from file, use comma as a delimiter
            vParsingVector.push_back(s); // Place value on vector
        }
        fData.close(); //Close file
    }

    else //Send error message if file cannot be opened and terminate program.
    {
        cout<< "Sorry Data (text) file could not be opened, check directory path, or the existence of the file and try again" << endl;
    }

    *iIterator = vParsingVector.size()-11; //Use iterator to declare velocity vector

    for( int i = 9; i<*(iIterator); i = i+10) // Cast string vector index to a double vector index, loops through the size of the vector.
    {
        double dbValue;
        stringstream ss;
        ss << vParsingVector[i];
        ss >> dbValue;



        vDataVector->push_back(dbValue);


    }

    *iIterator = vDataVector->size();
}



//Function to print calculated values into excel

void vfPrintValuesToExcel(int* iCorners, vector <double> *vCornerVelocityLosses, vector <double> *CornerDistances, vector <double> *vEnergyLosses,
                          vector <double> *vdHeatGeneratedByCorner, vector <double> *vdTimeonCorner,   vector <double> *vCornerInitialVelocity,
                          vector <double> *vCornerExitVelocity, vector <double> *vdCornerTempRise)
{
    ofstream ofAnswerFile; //create variable to open file

    ofAnswerFile.open("Corner_Calculations.csv"); //Print all values to this path

    if (ofAnswerFile.is_open())

    ofAnswerFile.clear(); //Clear contents before writing into file

    ofAnswerFile << " Corner , Entering Velocity (km/h), Exit Velocity (km/h), Velocity Loss (km/h), Length (m), Energy Loss (J)," <<
                 " Heat Generated (J), Time on Corner(s), Temperature Rise (C) " << endl; // makeheadings

    for (int i =0; i<*(iCorners); i++) // Print all values
    {
        ofAnswerFile << (i+1) << ", " << (vCornerInitialVelocity->at(i))*3.6 << ", " << (vCornerExitVelocity->at(i)*3.6)<< ", " << (vCornerVelocityLosses->at(i)*3.6) << ", " << CornerDistances->at(i) << ", "
                     << vEnergyLosses->at(i) <<  ", " << vdHeatGeneratedByCorner->at(i) << ", " << vdTimeonCorner->at(i) <<", " << vdCornerTempRise->at(i) << endl;
    }

    ofAnswerFile.close(); //close file
    cout << " ALL DATA PRINTED SUCCESSFULLY TO:" << endl << "Corner_Calculations.csv " << endl;
}

//Calculating current angular acceleration
double dbfCalcAngularAcceleration ( double *dbGearRatio, int iCurrentIndex, vector <double> *vdRPM)
{
         double dbCurrent_Motor_RPM = vdRPM->at(iCurrentIndex); // Need RPM and convert it for calculations

        dbCurrent_Motor_RPM = dbCurrent_Motor_RPM/60; //Get Revolutions per second
        double dbCurrent_Angular_Speed = dbCurrent_Motor_RPM*(2*M_PI); // Get angular speed per second by multiplying RPS per 2PI
        //Apply Gear Ratio
        dbCurrent_Angular_Speed = dbCurrent_Angular_Speed/(*dbGearRatio); // Get actual angular speed in front wheel by dividing by gear ratio

        return dbCurrent_Angular_Speed; //return value
}

double dbfCalcForce (double *dbGearRatio, int iCurrentIndex, vector <double> *vdTorque, double *dbFrontWheelRadius)
{
     double dbCurrent_Torque = vdTorque->at(iCurrentIndex); // Get initial torque of specific corner
     if (dbCurrent_Torque<0)
     {
         dbCurrent_Torque = dbCurrent_Torque*-1;
     }

        dbCurrent_Torque = dbCurrent_Torque/(*dbGearRatio); // Pass Torque through gear ratio to get torque at the front wheel
        dbCurrent_Torque = dbCurrent_Torque/(*dbFrontWheelRadius); // Divide by front wheel radius to calculate force

        return dbCurrent_Torque; //return value
}

double dbfCalculateForce (int z, vector <double> *vdVelocity, double* bike_mass, vector <double> *vdDistance, double* ke)
{
    int prev = z-1;
    double step_energy_loss = ((vdVelocity->at(prev)*vdVelocity->at(prev))-(vdVelocity->at(z)*vdVelocity->at(z)))*(*bike_mass)*0.5;
    *ke = 0.5*step_energy_loss;
    double step_distance = (vdDistance->at(z)-vdDistance->at(z-1));
    double step_force = step_energy_loss/step_distance;

   //x cout<< "Prev Velocity = " << vdVelocity->at(prev) << " Current Velocity: "<< vdVelocity->at(z) << " step ke: " << step_energy_loss << endl;
    return step_force;

}

//Function to calculate Step Temperature Rise for the brake disc
double dbfStepHeatConductedOut(double* dbBrakeDiscK, double* dbBrakeDiscThickness, double* dbPadArea, double* dbHeatGenerated)
{
    double HeatOut = ((*dbHeatGenerated)*(*dbBrakeDiscThickness))/ ((*dbBrakeDiscK)*(*dbPadArea)); //Apply calculation
    return HeatOut;
}

void vfProcessBruntTestData(vector <double> *Time, vector <double> *BrakeTemp, vector <double> *Speed, string *strFilePath)
{
    fstream fTestData; //fstream variable to open variable

    fTestData.open(strFilePath->c_str());
    if (fTestData.good())
    {
        while (!fTestData.eof())
        {
            string sTime; //for direct input
            string sTemp;
            string sSpeed;
            string s; // for parsing un-needed values
            string si;
            stringstream ssTime; //to convert values from strings to doubles
            stringstream ssTemp;
            stringstream ssSpeed;

            double dbTime; //double values to place in vector
            double dbTemp;
            double dbSpeed;


            getline(fTestData, si, '\n'); // need to read newline character
            stringstream ssLine(si);
            getline(ssLine, sTime, ',');
            getline(ssLine, sTemp, ',');
            getline(ssLine, sSpeed, ',');

            ssTime << sTime; //convert string to double
            ssTime >> dbTime;
            Time->push_back(dbTime);// place it on required rate

            ssTemp << sTemp;
            ssTemp >> dbTemp;
            BrakeTemp->push_back(dbTemp);

            ssSpeed << sSpeed;
            ssSpeed >> dbSpeed;
            Speed->push_back(dbSpeed);


        }
        fTestData.close();
    }

    else
    {
        cout << "Sorry unable to properly process data, check file path and try again \n";
    }
}

void vfPrintTestData(int *iSize, vector <double> *vTime, vector <double> *vTemp, vector <double> *vSpeed, vector <double> *vBrakeTemp)
{
    ofstream ofAnswerFile; //create variable to open file

    ofAnswerFile.open("INITIAL_GUESSES.csv"); //Print all values to this path

    if(ofAnswerFile.is_open())
    {

    ofAnswerFile.clear(); //Clear contents before writing into file

    ofAnswerFile << " Time (s) ,Temperature(C), Velocity (km/h), Estimated Temperature(C)  " << endl; // makeheadings

    for (int i =0; i<*(iSize); i++) // Print all values
    {
        ofAnswerFile << vTime->at(i) << ", " << vTemp->at(i) << ", " << vSpeed->at(i)<< "," << vBrakeTemp->at(i) << endl;
    }

    ofAnswerFile.close(); //close file
    }

    else
    {
        cout << "Unable to produce estimated plot data" << endl;
    }

}

void vfPrintFinalTempData(int iSize, vector <double> *vDistance, vector <double> *vTemp, vector <int> *vVelocity)
{
    //int iSize = vDistance->size();
    ofstream ofAnswerFile; //create variable to open file

    ofAnswerFile.open("TT_Estimate.csv"); //Print all values to this path

    if (ofAnswerFile.is_open())
    {

    ofAnswerFile.clear(); //Clear contents before writing into file

    ofAnswerFile << " Distance (m) ,Temperature(C), Velocity (km/h)" << endl; // makeheadings

    for (int i =0; i<(iSize); i++) // Print all values
    {
        ofAnswerFile << vDistance->at(i) << ", " << vTemp->at(i) << ", " << vVelocity->at(i) << endl;
    }

    ofAnswerFile.close(); //close file
    }

    else
    {
        cout << "Unable to create final temperature file" << endl;
    }

}
