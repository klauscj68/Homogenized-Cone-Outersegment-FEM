This directory contains all the files needed to run the homogenized cone model.

Step 1:
A user first goes to the common/ folder and opens data.m . This file contains all
the geometric, biochemical, and numerical parameters used throughout the simulation.
Similarly, the number of photons deposited in the COS and where are also specified
by this file. 

Step 2:
After editing the data.m file, a user proceeds to the id/ folder. Run main_id.m .
This program simulates how detection of the photons at a special disc generates
G-protein and Phosphodiesterase, E_st. This data is saved in the file E_st.mat .

Step 3:
Once E_st.mat is generated, a user proceeds to the cyto/ folder.  Now run main_cyto.m .
This program now simulates how E_st on a disc produces volumic changes in cGMP and 
Ca2+ throughout the cytosol. When the program concludes, all the relevant data is 
saved in cyto.mat .  

Note: 
Solutions of cGMP and Ca2+ are stored as coefficients on a nodel basis.  
The ordering of the nodal basis is described in the file cyto/gdl.m .

Note:
Files ending .mex* in the elements/ folder were compiled from matlab from C code.  It is 
possible that your computer may need to recompile the mex code from the C files
to run properly.
