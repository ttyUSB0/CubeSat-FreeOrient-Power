# Estimation of the CubeSat's available energy for free-orientation scenario

The power budget of CubeSat is difficult to estimate, especially when  the law of its orientation while moving in orbit is unknown. If the orientation system is absent, failed, or an orientation failure has occurred, then the satellite usually begins to rotate; in this case, the axis of its own rotation can be any. The amount of available energy in free-orientation mode determines how quickly it can restore its functionality and whether it can at all. Thus, the available energy can be obtained using a statistical estimate. 

This article (and code!) provides estimates for popular CubeSat designs. To estimate the energy, one need to substitute the parameters of the solar cells used (at maximum power point), as well as the fill factor of the solar array panels.

For more details see article in pdf file.
If you found an error in math or physics, please email me to a.t.lelekov АТ  yandex ДОТ ru.

# Project structure
Code divided by two parts. First calculates statistical parameters for cubesat 1U; second apply weights for satellite sides and visualizing results.
I execute code (python3) in IPython console in Spyder (Anaconda) in Debian10. These are scripts, and must be started in IPython console in hand mode. Sorry, I haven't prepare code for Windows platform, check file operations.

## Statistics.py
I provide results of statistics calculation in file "cube.pickle", but you can get it yourself with desired parameters:
    a = np.linspace(0, 2 * np.pi, 180) defines tolerance of integral. Change 180 to 360 to get more precise value.

    NPoints = 250 (int) defines number of vectors v on sphere (see article). You can increase it to get result more precise.

    db = 10 (int) defines delta in rotation angle, when cube rotates over v from origin to initial position. You can lower it to get result more precise.

To run script you must provide path to project folder (change variable "folder", sorry :). 
Then script will make all math, descripted in article. 

It will generate a set of vectors v, uniformly distributed in all possible directions. Then it will rotate cube (represented by its normals) over all vectors and integrate illumination of each normals over this turn -- in function NormalsCos. It returns six values of integral, for each normal.
Finally, script saves resulting numpy array of integral to file "cube.pickle". To get results, go to Visualize.

## Visualize.py
Firstly check variable "folder". Then execute block "Read data", it will load results.
Next, execute code with definition of Sats - for 1U/2U/3U structure, which define weights for normals.
Finally, execute remaining code to plot dfs and get statistical values in cout.


You can use this code under MIT license, 
with copyright (c) 2020 A.T.Lelekov