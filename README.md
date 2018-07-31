# UoN Racing: IoM TT Zero Superbike Brake Validation
***

This repository contains the now completed code and documentation for the complete brake validation, which was part of our project for the University of Nottingham's Racing Team Superbike that achieved second place in the TT Zero Division of the prestigious race. 

## Code and Specifications
***
All the available code is in the ```TT_Zero_Validation_Code```. It was coded in C++ using CodeBlocks (an [open source IDE](http://www.codeblocks.org/)) using the GNU GCC Compiler (also open source!). The project has around ~700 lines of code. There are attached documents which explain the formulae  and the methodology used to complete the brake validation, and the main algorithm. 

## Why validate the brakes?
***
Validating the brake discs in any vehicle is very important, in a racing vehicle it is paramount. Electric racing bikes are a lot heavier than standard (petrol) racing bikes, in fact UoN's superbike weighs ~275kg (about 100kg more than the average petrol racing bike) and go to speeds up to 280 km/h. This causes some issues from a safety and engineering standpoint, as superbikes being 2-wheel vehicles only has 2 brake discs (on the front wheel) to dissipate all the heat caused by braking (cars have 4). This means a huge amount of heat needs to be dissipated in comparison to other vehicles (even other racing bikes), and thus finding the right brake configuration becomes crucial to prevent crashes which can lead to fatalities. 

The program developed takes in a lot of parameters (e.g. brake disc size, thickness, material, heat transfer coefficient, thermal diffusivity, cooling rates) and uses heating and cooling equations to monitor the brake disc temperature, making sure it does not go above its maximum service temperature to keep the rider safe and maximise racing performance. The final results can be seen in the spreadsheet in the repository called ```FINAL RESULTS```. General characteristics of the brake forces, and other attributes can be found in the spreadsheet in the repository called ```Full Brake Calculated Data```. Explanation in detail of the calculations used is found on the pdf ```Formulaes and Calculations Explanation``` and the final graph by itself can be seen in ```Brake Performance```

## Solving the problem
*** 
The basis of the solution is raw data recorded by a datalogger attached to the superbike, when it competed last year. The datalogger recorded at a sampling rate of about 80 times a second and provided thousands of values of useful data (i.e. distance, velocity, motor rpm). 

An assumption was made that in racing conditions the bike would only brake when in corners. Hence, the first portion of the algorithm sorts through the thousands of values for velocity, locates continuous loss of velocity and tags them as corners. The corner classifications are then sanitised (dummy values such as velocity losses of <1 m/s are removed), and using kinematic energy equations the energy needed to be dissipated by the brake discs per corner was calculated. 

From the brake disc values, and the bike's total mass, the brake disc temperature rise per corner was calculated. While very useful in brake validation, did not give the full picture if the brake disc was suitable for the race. 

Therefore, the next step was to calculate the cooling the brake disc would experience in between corners, this was calculated by getting cooling coefficients by experimentation, (attaching an infrared sensor to the superbike) and therefore integrating and applying Newton's Law of Cooling, which allowed us to know the estimated brake disc temperature at any given time. By knowing the disc brake temperature at any time, we can validate the brakes by checking if at any given point(s) it exceeds its lowest maximum service temperature. Hence fully validating the brakes.




