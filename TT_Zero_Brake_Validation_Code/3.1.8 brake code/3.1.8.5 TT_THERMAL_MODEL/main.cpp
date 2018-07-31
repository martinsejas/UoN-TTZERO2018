#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "thermal_model_functions.h"
#include "calculate_thermal_coefficients.h"

// Import file input and reading header files
using namespace std;


int main()
{

    vector <double> vBrakingCoefficients;
    vector <double> vCoolingCoefficients;

    vfCalculateThermalCoefficients(&vBrakingCoefficients, &vCoolingCoefficients);





    /* Explain to user how it works*/

    int iVectorIteratorV; //Integer variable to iterate through velocity vector
    int iVectorIteratorD; //Same  but for Distance vector
    int iVectorIteratorT; // Torque
    int iVectorIteratorRPM; //RPM


    string velocity_file_path = "VELOCITY.txt"; //Get the required file paths
    string distance_file_path = "DISTANCE.txt";
    //string torque_file_path = "C:\\Users\\MartinHugo\\Desktop\\Bike_Sim\\data_export\\Torque_Plot_Data_NoRB.txt";
    //string RPM_file_path = "C:\\Users\\MartinHugo\\Desktop\\Bike_Sim\\data_export\\RPM_Plot_Data_NoRB.txt";

    vector <double> vdVelocity; // Assign a vector to each variable needed
    vector <double> vdDistance;
    vector <double> vdTorque;
    vector <double> vdRPM;
    vector <double> vdCornerDistances;
    vector <double> vdAngularSpeedByCorner;
    vector <double> vdInitialTorqueByCorner;
    vector <double> vdInitialForcebyCorner;

    vfReadTextIntoVector(&velocity_file_path, &vdVelocity, &iVectorIteratorV); //Call function to set text (data) file into a vector variable
    vfReadTextIntoVector(&distance_file_path, &vdDistance, &iVectorIteratorD);
//    vfReadTextIntoVector(&torque_file_path, &vdTorque, &iVectorIteratorT);
    //vfReadTextIntoVector(&RPM_file_path, &vdRPM, &iVectorIteratorRPM);

    // Time to check for corners, count them and declare velocity loss (energy loss) per corner as well as corner length.
    cout << "I GOT HERE (PAST VECTOR ASSIGNMENTS) " << endl;



    vector <string> vCornerIndexes; //Declare a vector to hold the start and end of corners indexes ( to find distances later)
    vector <double> vCornerEnergyLosses;//Vector holding kinetic energy losses for each corner
    vector <double> vCornerEnteringVelocity; // vector holding initial entering velocities
    vector <double> vCornerExitVelocity; // vector holding exiting velocities
    vector <double> vCornerVelocityLosses; // Vector holding velocity losses for each corner
    vector <double> vdTimeOnCorner; // vector holding time in each corner
    vector <double> vdPressureOnCorner; // Declare vector to hold pressure caused by brake callipers on brake disc by corner
    vector <double> vdHeatGeneratedByCorner;// vector holding heat generation per corner
    vector <double> vdCornerTempRise;
    int iStartIndex; // integer variable to hold the start of a corner
    int iEndIndex; // integer variable to hold the end of a corner
    double dPrevious = 0; // integer variable to compare previous value to current value to determine if a corner has started
    double dCurrent = 0; // integer variable that iterates through vector
    double dStartCorner; // double variable to hold the velocity at the start of the corner
    double dEndCorner; // double variable to hold the velocity at the end of the corner
    bool bInCorner = false; //boolean to determine if current index is inside a corner or not
    double total_bike_mass = 380; // From simulation, 380 total bike mass
    double dbTempRiseFactor = 3.07; // Correction factor
    double dInitialTemp = 20.0;

    vector <int> vCornerStart;
    vector <int> vCornerEnd;

    cout << "IVECTORITERATOR: " << iVectorIteratorD << endl;


    for(int i = 1; i<iVectorIteratorD; i++) // loop through entire array, start with i = 1 because there is the dPrevious variable
    {
        dPrevious = vdVelocity[i-1]; //Assign values
        dCurrent = vdVelocity[i];

        if(!bInCorner && dCurrent<dPrevious) // If index is not in a corner check to see if there is a velocity drop
        {
            dStartCorner = dPrevious; // if there is, it means the start of a corner, save velocity value before braking started
            iStartIndex = i-1; // Save the index of the start of the corner
            bInCorner = true; // we are now in a corner
        }

        if (bInCorner && dCurrent>dPrevious) // If index is in a corner but velocity has increased, the driver has exited the corner
        {
            dEndCorner = dPrevious; //Save the velocity at the end of the corner
            iEndIndex = i-1; //Save index indicating end of corner
            stringstream sstream; // Declare sstring to convert integer values into a string (for parsing later)
            sstream << iStartIndex << "," << iEndIndex; // save indexes using a comma as a delimiter
            string s = sstream.str(); // set string value to indexes
            vCornerIndexes.push_back(s); // place in vector holding indexes
            vCornerVelocityLosses.push_back(dStartCorner-dEndCorner);

            double energy_loss = ((vdVelocity[iStartIndex]*vdVelocity[iStartIndex])-(vdVelocity[iEndIndex]*vdVelocity[iEndIndex]))*(total_bike_mass/2); //place calculated energy loss
            vCornerEnergyLosses.push_back(energy_loss); // Place velocity loss in velocity vector
            bInCorner = false; // we are now out of a corner
        }
    }


    cout << " GOT SAFELY OUT OF CORNER CALCS!!!" << endl;
    cout << " THERE THIS: " << vCornerIndexes.size() << "MANY CORNERS" << endl;



    for (unsigned int i = 0; i<vCornerVelocityLosses.size(); i++ ) //iterate through velocity vector and erase phony values (m/s < 1)
    {
        if (vCornerVelocityLosses[i] < 1)
        {
            vCornerVelocityLosses.erase(vCornerVelocityLosses.begin()+i);
            vCornerEnergyLosses.erase(vCornerEnergyLosses.begin()+i);
            vCornerIndexes.erase(vCornerIndexes.begin()+i); // also erase from index vector
        }
    }


    cout << "GOT OUT OF PHONY VALUES!!" << endl;
// calculating total corner distance

    int iCorners = vCornerIndexes.size();


    for (int i = 0; i<iCorners; i++) //Parse through corners and make a vector of the distances of each corner
    {
        int iStartDistanceIndex;// Variable to index start of corner
        int iEndDistanceIndex; // End of Corner

        stringstream ss(vCornerIndexes[i]); // stringstream to convert string to int from corner indexes


        ss >>iStartDistanceIndex; // String to int (first value of index)
        ss.ignore(); //Ignore delimiter
        ss>> iEndDistanceIndex;// String to (final value of index)

        vCornerStart.push_back(iStartDistanceIndex);
        vCornerEnd.push_back(iEndDistanceIndex);

        cout << "vCornerStart[" << i << "] is: " << vCornerStart[i] << endl;
        cout << "vCornerEnd[" << i << "] is: " << vCornerEnd[i] << endl;

        double corner_distance = (vdDistance[iEndDistanceIndex]-vdDistance[iStartDistanceIndex]); // Calculate specific corner distance
        vdCornerDistances.push_back(corner_distance); //Place it on a respective vector

    }

    //cout << "Total Corner Distances = " << vdCornerDistances.size() << endl; //Output distances



    vector <double> vdDeAccelerationRatePerCorner; //DeaccelerationRaterPerCorner

    for( int i =  0; i<iCorners; i++)
    {
        double deAcceleration_Rate = vdCornerDistances[i]/vCornerVelocityLosses[i];
        vdDeAccelerationRatePerCorner.push_back(deAcceleration_Rate);
    }


    double dbGearRatio = 83/17; //Obtained from workshop measurements/calculations
    double dbFrontWheelRadius = 0.6; // metres, measured in workshop
    double dbTotal_Piston_Area = 0.007636; // 0.1885 Calculated piston area by measurements in workshop
    double datalogger_unit_time = 0.01112; //Time in seconds between each index interval
    double dbEffective_Radius_Area = 0.0355739;
    double dbEffective_Radius = 0.106412;
    double dbCoefficientFriction = 0.41;
    double dbHeatPartition = 0.88;
    double dbBrakeDiscK = 40; // thermal conductivity in for ductile iron 32.3 Wm-1K-1
    double dbBrakeDiscThickness = 0.00725; // 0.00725 thickness in metres (delta x)
    double dbPadArea = 0.004011; //In metres squared, ESTIMATE
    double dbHeatCapacity = 480; // 506 for ductile iron
    double dbBrakeDiscArea = 0.02673078; // m^2
    double dbBrakeVolume = 0.0001978078; // m^3
    double dbBrakeDiscMass = 1.4; // in Kg
    double dbBrakeDiscDensity = 7150;// kg/m^3
    double dbThermalDiffusivity = dbBrakeDiscK/(dbBrakeDiscDensity*dbHeatCapacity);//0.0000957 1.2*10^-5


    // double dbBrakeSurfaceArea =



    //vdRPM, get RPM for start of corners, divide to get angular speed, then apply gear ratio,
    //and then use initial RPM and with average velocity per corner calculate braking angle



    //torque needed to calculate radial pressure, Torque = Fd, torque on wheel calculated by gear ratio, wheel diameter known.
    // Sliding force can be modelled equal to direct force due to Newton's second law.
    // Divide force over the area of the 6 pistons to get pressure
    // Calculate heat on pad

    int calc_brake_indexes = 0;
    for (int i = 0; i<iCorners; i++) //Parse through corners match RPM  and Torque values for each
    {
        int iSIndex;//starting index
        int iEIndex;//ending index

        stringstream ss(vCornerIndexes[i]);


        ss >>iSIndex; // Match initial index value for specific corner
        ss.ignore();
        ss>> iEIndex;
        calc_brake_indexes = calc_brake_indexes+(iEIndex - iSIndex);
        double total_heat_generated = 0;
        double total_kinetic_energy = 0;
        double kinetic_energy = 0;
        double step_heat = 0;
        double internal_heat_in = 0;
        double step_temp = 0;
        double total_temp = 0;
        double total_internal = 0;
        double total_heat_out = 0;
        double time_on_corner = (iEIndex-iSIndex)*datalogger_unit_time;

        for( int z = iSIndex; z<=iEIndex; z++) // perform for loop inside loop to simulate integration
        {

            double current_energy= dbfCalculateForce(z,&vdVelocity, &total_bike_mass,&vdDistance, &kinetic_energy);
            total_kinetic_energy = total_kinetic_energy+kinetic_energy;


        }
        total_heat_generated = total_kinetic_energy;

        total_temp = total_kinetic_energy/(dbBrakeDiscMass*dbHeatCapacity);
        vdCornerTempRise.push_back(total_temp/3.7);

        //double total_temp_rise = (total_heat_generated-total_heat_out)/(dbBrakeDiscMass*dbHeatCapacity);
        double total_temp_rise = (0.53*total_heat_generated)*pow(((time_on_corner)/(dbBrakeDiscK*dbBrakeDiscDensity*dbHeatCapacity)), 0.5);
        //cout << " total heat generated is: " << total_heat_generated << " time on corner is: " << time_on_corner << " pow term is: "<< pow(((time_on_corner)/(dbBrakeDiscK*dbBrakeDiscDensity*dbHeatCapacity)), 0.5)<< endl;
        vCornerEnteringVelocity.push_back(vdVelocity[iSIndex]);
        vCornerExitVelocity.push_back(vdVelocity[iEIndex]);
        //cout << "Corner: " << (i+1) << " I Vel: " << vdVelocity[iSIndex] << " F Vel: " << vdVelocity[iEIndex] << " Temp_UF: " << total_temp_rise << endl << endl;
        vdHeatGeneratedByCorner.push_back(total_heat_generated); // total heat generated calculated
        vdTimeOnCorner.push_back((iEIndex-iSIndex)*datalogger_unit_time); //calculating time on corner
    }

    //==============MAIN LOOP========================
    //==============IGNORE VALUE UNTIL FIRST CORNER===============
    //==============CHECK IF NEXT CORNER AT EVERY ITERATION=============
    //==============IF CORNER ADD TEMP RISE AND SHIFT INDEXES===========
    //==============STORE GUESSED TEMP IN A NEW VECTOR =================
    //==============VECTOR WILL HAVE THOUSANDS OF VALUE PRINT EVERY 10 OF THEM=========
    double dEstimatedTemp = dInitialTemp;
    double dPrevTemp = dEstimatedTemp;
    vector <double> vEstimatedTemps; //to store values
    vector <int> viVelocity;
    int iC = 0; //iterator for corners
    int piC = iC;
    int total_brake_indexes;
    int brake_count = 0;
    double dbCoolingCoefficient = 0.00413;
    double unit_time = 10*datalogger_unit_time;

     vector <double> dvHvalues;

    for (int i = 0; i<300; i++)
    {
        double v = i/3.6;
        dvHvalues.push_back(((dbCoolingCoefficient*dbBrakeDiscMass*dbHeatCapacity)/dbBrakeDiscArea) + pow(v,0.6));
        cout << dvHvalues[i] << endl;
    }

    for (unsigned int i =0; i<vdVelocity.size(); i++) // convert to km/h
    {
        int vi = (int)((vdVelocity[i]*3.6) + 0.5); // got to truncate it
        viVelocity.push_back(vi);
    }

    for (unsigned int i = 0 ; i< viVelocity.size()-1; i++)
    {

        if (i < vCornerStart[0])// Ignoring values until first corner
        {

            dEstimatedTemp = dInitialTemp;
            vEstimatedTemps.push_back(dEstimatedTemp);
            dPrevTemp = dEstimatedTemp;
        }

        else if (i == vCornerStart[iC]) //checking if we are in a corner
        {
            brake_count++;
            //if we are in a corner the temperature has to be added, and needs to be split evenly
            int rang = vCornerEnd[iC]-vCornerStart[iC];
            double Temp_Increment = vdCornerTempRise[iC]/( (double) rang);
            //if (i >3133 && i< 3136) cout << " Total TEMP CORNER RISE: " << vdCornerTempRise[iC]  << endl;
            rang = rang+i;
            for (; i < (rang); i++)
            {
                dEstimatedTemp = dEstimatedTemp+Temp_Increment; //incremental rise
                vEstimatedTemps.push_back(dEstimatedTemp);
                dPrevTemp = dEstimatedTemp;
            }
            piC = iC;
            iC++;
        }


        else if (i != vCornerStart[iC] && i>vCornerStart[0])
        {

           /* int diff = viVelocity[i]-viVelocity[(i+1)];
            diff = abs(diff); // get velocity diff
            double dTdiff = dEstimatedTemp - dInitialTemp; */

            double time_passed = (i - vCornerEnd[piC])*unit_time;
            double heat_coefficient = dvHvalues[viVelocity[i]];
            double e = 2.71828;
            double expo = -(heat_coefficient*dbBrakeDiscArea/(dbBrakeDiscMass*(dbHeatCapacity/10))*time_passed);

            dEstimatedTemp = dInitialTemp + pow(e,expo)*(dPrevTemp-dInitialTemp);

           /* if (i < 447 && i> 439)
            {
                cout << "Velocity diff: " << diff << " Estimated Temp: " << dEstimatedTemp << " dInitial Temp: " << dInitialTemp << " dTdiff: " << dTdiff << "Coefficient: " << vCoolingCoefficients[viVelocity[(i+1)]] << endl;

            }

                dEstimatedTemp = (vCoolingCoefficients[ viVelocity[(i)] ]*dTdiff); // add the coefficients to the temp
                */
                if (dEstimatedTemp < dInitialTemp)
                {
                    dEstimatedTemp = dInitialTemp+0.1;
                }


                vEstimatedTemps.push_back(dEstimatedTemp);


        }

    }



    unsigned int temp_size = vEstimatedTemps.size();
    cout << "brake corner count: " << brake_count << endl;
    cout << " corner array size" << iCorners << endl;
    cout << "vdVelocity size is: " << vdVelocity.size() << endl;
    cout << "viVelocity size is: " << viVelocity.size() << endl;
    cout << "Temp size is: " << vEstimatedTemps.size() << endl;
    cout << " DISTANCE V IS: " << vdDistance.size() << endl;

    vfPrintFinalTempData(temp_size, &vdDistance, &vEstimatedTemps, &viVelocity);

    // vfPrintFinalTempData(int iSize, vector <double> *vDistance, vector <double> *vTemp, vector <int> *vVelocity);

    vfPrintValuesToExcel(&iCorners,&vCornerVelocityLosses,&vdCornerDistances,
                         &vCornerEnergyLosses,&vdHeatGeneratedByCorner,&vdTimeOnCorner,&vCornerEnteringVelocity, &vCornerExitVelocity, &vdCornerTempRise);

    double total_time_on_corners = 0;
    double total_corner_length = 0;
    for ( unsigned int i = 0; i< vdTimeOnCorner.size(); i++)
    {
        total_time_on_corners = total_time_on_corners+vdTimeOnCorner[i];
        total_corner_length = total_corner_length+vdCornerDistances[i];

    }


    vector <double> vdNewTimeOnCorner;
    double total_new_time = 0;
    for (unsigned int i = 0; i<vdCornerDistances.size(); i++)
    {
        double dbAverageVelocity = vdCornerDistances[i]/vdTimeOnCorner[i];
        dbAverageVelocity = dbAverageVelocity+0.75;
        vdNewTimeOnCorner.push_back((vdCornerDistances[i]/dbAverageVelocity));
        total_new_time = total_new_time+vdNewTimeOnCorner[i];
    }

    return 0;
}
