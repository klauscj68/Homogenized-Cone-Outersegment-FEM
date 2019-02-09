This folder contains all the files needed to run the nonhomogenized cone model.

Step 1:
A user opens the data.m file and adjusts all the relevant geometric, biochemical,
and numerical parameters here as well as specifying where photons are to be 
deposited. 

Step 2: 
A user should decide whether the simulation needs to also build the mesh and 
mass and stiffness matrices or if instead these quantities may be preloaded
because they were built before. If they need to be built, one opens master.m
and sets flag_build_mesh = true;. Next one opens assembla.m and sets 
flag_assemble = true;. If instead they can be preloaded, set these flags to false
and make sure that the mesh and matrix .mat files haven't had their names changed.

Step 3:
One runs master.m . This will generate the cyto.mat and cyto_sol.mat where all the 
data is stored.  cyto.mat stores all the output variables except the solution of 
cGMP and Ca2+ stored over the 3d mesh.  These are instead saved in the separate 
cyto_sol.mat file owing to their size.

Note:
Files ending .mex* were compiled from matlab from C code.  It is 
possible that your computer may need to recompile the mex code from the C files
to run properly.
