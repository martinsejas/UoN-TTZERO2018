#include <iostream>
#include <vector>
#include <string>
#include "thermal_model_functions.h"

using namespace std;

void vfCalculateThermalCoefficients(vector <double> *vBrakingCoefficients, vector <double> *vCoolingCoefficients)
{
    vector <double> IR_Time; // Set vectors for retrieving data
    vector <double> IR_BrakeTemp;
    vector <double> IR_Velocity;


    vector<double> vRaceDiscTemp;

//============NOTE: HOW THE COEFFICIENTS WORK===========================================================================
// IF WE WANT THE TEMPERATURE DIFFERENCE WHILE BRAKING (TEMP RISE AS DISC GETS HOTTER) FROM SAY 180km/h to 179 km/h
// THE LOWEST VALUE IS THE INDEX FOR THE VECTORS AKA lets say the TEMP RISE from 180km/h to 179 km/h be 0.9 C
// YOU WOULD REFERENCE vBrakingCoefficients[179] and add it to the Brake Disc TEMP
// SAME THING FOR COOLING COEFFICIENTS!!
// =========== END NOTE ================================================================================================

    string IR_Path = "BRUNTI_TEST_DATA.csv";


    vfProcessBruntTestData(&IR_Time, &IR_BrakeTemp, &IR_Velocity, &IR_Path); //Function that will read .csv file accurately
    int iSize = IR_Time.size(); // Get the total value of the arrays



    for (int i = 0; i<300; i++ ) // initialize vectors and array
    {

        vBrakingCoefficients->push_back(0);
        vCoolingCoefficients->push_back(0);
    }



    double dInitialTemp = 16.0;
    int next_i; // next index
    int n; // variable used to get velocity difference Zv - Z+1v = n
    double pT; //previous temp
    double dT; // temperature difference
    double dTroom; // difference between ambient to origintal
    double dStack_T; // for when there is the same vel
    int iCurrentVelocity; // for calcs
    int iNextVelocity; // for calcs
    bool bBraking = true; // to determine if cooling or not ! POSSIBLY CHANGE TO bMOVING

for (int i = 0; i<1000; i++)
{


    for (int z = 0; z < iSize-1; z++) // main array to get coefficients, starting at 257 because thats when first brake starts
    {
        next_i = z+1;
        iCurrentVelocity = IR_Velocity[z]; // getting velocities
        iNextVelocity = IR_Velocity[next_i];

        n = iCurrentVelocity - iNextVelocity; // see velocity difference, three cases apply here

        if (n>0) // If n is positive it means it is currently undergoing breaking
        {
            bBraking = true; //set boolean for positive as we are indeed braking
            dT = IR_BrakeTemp[next_i]-IR_BrakeTemp[z]; // get temperature differences from indexes
            dTroom = IR_BrakeTemp[next_i]-dInitialTemp; // get temperature difference from ambient
            //stack conditional
            if(dT<0) //if bogus values
            {
                dT = 0;
                pT = dT;
            }
            else
            {
                dT = dT+dStack_T; // add stack if needed
                pT = dT;
                dT = dT/( (double)n); // split if equally if there is a spread of velocities

                dStack_T = 0;
            }

    // CONSIDER FOR RaNGED VALUES
            for (; iNextVelocity<iCurrentVelocity; iNextVelocity++) // assign temperature differences
            {
                double y = (dT+IR_BrakeTemp[next_i])/dTroom;
                if(vBrakingCoefficients->at(iNextVelocity) != 0) // if not empty get the average
                {
                    vBrakingCoefficients->at(iNextVelocity) = (vBrakingCoefficients->at(iNextVelocity)+y)/2.0;
                }
                else
                {

                    vBrakingCoefficients->at(iNextVelocity) = y;
                }

            }


        }
        else if (n<0) // if n is negative it means we are accelerating
        {
            bBraking = false; //not braking
            dT = IR_BrakeTemp[next_i]-IR_BrakeTemp[z]; // get temperature differences from indexes
             dTroom = IR_BrakeTemp[next_i]-dInitialTemp;
            if(dT>0) //if bogus values
            {
                dT = 0;
                pT = dT;
            }
            else
            {
                dT = dT+dStack_T; // add stack if needed
                pT = dT;
                dT = dT/( (double)n); // split if equally if there is a spread of velocities
                dStack_T = 0;
            }

            for (; iCurrentVelocity<iNextVelocity; iCurrentVelocity++)
            {
                double y = (dT+IR_BrakeTemp[next_i])/dTroom;
                if (vCoolingCoefficients->at(iCurrentVelocity) != 0)
                {

                    vCoolingCoefficients->at(iCurrentVelocity) = (vCoolingCoefficients->at(iCurrentVelocity)+y)/2.0;
                }
                else
                {

                    vCoolingCoefficients->at(iCurrentVelocity) = y;
                }
            }

        }

        else if (n == 0)
        {
            if (iNextVelocity == 0)
            {
                dStack_T = 0;
            }

            else if(pT != 0)
            {
                dStack_T = IR_BrakeTemp[next_i]-IR_BrakeTemp[z];
                if (bBraking == true)
                {
                    dStack_T = dStack_T*-1;
                }

            }
            else
            {
                dStack_T = dStack_T;
            }

        }
        else
        {
            cout << " ERROR: UNABLE TO SORT COEFFICIENTS" << endl;
        }

    }

     for (unsigned int it = 0; it<300 ; it++)
    {

        if (vCoolingCoefficients->at(it) < 1 || vCoolingCoefficients->at(it) > 1.2)
        {
            if (it <80) vCoolingCoefficients->at(it) = 1.03;
            if (it <= 150 && it>80) vCoolingCoefficients->at(it) = 1.05;
            if (it> 150) vCoolingCoefficients->at(it) = 1.08;
        }

            //cout << "Cooling Coefficient[" << it << "] is: " << vCoolingCoefficients->at(it) << endl;

    }

}




    //NOW YOU APPLY IT TO THE TEST DATA
    // ALWAYS REFERENCE SMALLER DATA


    double dEstimatedTemperature = dInitialTemp;
    for (unsigned int it = 0; it < IR_Velocity.size(); it++) // GET INITIAL ESTIMATES
    {
        int itt = it+1;
        int diff = IR_Velocity[it]-IR_Velocity[itt]; // get velocity diff
        dTroom = IR_BrakeTemp[it]-dInitialTemp;

        if (diff>0) // if we are braking
        {
            for (int x = 0; x< diff; x++) // add coefficients
            {
                if(dTroom>1)
                {
                  //dEstimatedTemperature = dEstimatedTemperature+(vBrakingCoefficients[IR_Velocity[it]+x]*dTroom); // add the coefficients to the temp
                  dEstimatedTemperature = (vBrakingCoefficients->at(IR_Velocity[it]+x)*dTroom);
                }

            }
            vRaceDiscTemp.push_back(dEstimatedTemperature);
        }

        else if (diff<0) // if we are accelerating
        {
            diff = diff*-1; // make the diff positive

            for (int x = 0; x< diff; x++) // add coefficients
            {
                if (dTroom>1)
                {

                    dEstimatedTemperature = (vCoolingCoefficients->at(IR_Velocity[it]+x)*dTroom); // add the coefficients to the temp
                }


            }

            vRaceDiscTemp.push_back(dEstimatedTemperature);
        }

        else
        {
            vRaceDiscTemp.push_back(dEstimatedTemperature);
        }
    }



    vfPrintTestData(&iSize,&IR_Time,&IR_BrakeTemp, &IR_Velocity, &vRaceDiscTemp);// print values

}
