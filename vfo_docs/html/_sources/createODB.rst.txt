.. include:: sub.txt

=========================
 Create Output Database 
=========================

.. function:: vfo.createODB(ModelName, <LoadCaseName>, <Nmodes=0>, <deltaT=0.0>, <recorders=[]>)

   This command creates an Output Database for the active model with an option to save a specific load case output.
   The command **must** be called while the model is built, but before the main analysis is run.
   An output database with name ModelName_ODB is created in the current directory with node and element information
   in it. See the example below.
   
   Input arguments are as follows,
   
   ==========================  ===============================================================================
   ``ModelName``     |str|      Name of the model the user wants to save database with.
   ``LoadCaseName`` |str|      Name of the subfolder to save load case output data.(Optional)
   ``Nmodes``        |int|      Number of modes to be saved for visualization.(Optional)
   ``deltaT``        |float|    Timesteps at which output to be saved. Default is 0.0. (optional)
   ==========================  ===============================================================================

Here is a simple example:

::

   vfo.createODB("TwoSpan_Bridge", Nmodes=3)

   
The above command will create,

* a folder named **TwoSpan_Bridge_ODB** containing the information on the nodes and elements to visualize the structure in future without using a OpenSeesPy model file.
* a sub-folder named **ModeShapes** containing information on modeshapes and modal periods.
* no output from any loadcase will be saved saved in this case even if the analysis is run in the script since there is no argument for "LoadCaseName" is provided.
   

Here is another example command to be used when recording data from a load case.

::

   vfo.createODB("TwoSpan_Bridge","Dynamic_GM1", Nmodes=3, deltaT=0.5)

   
The above command should be used right before running a load case analysis in the OpenSeesPy script and will create,

 * a folder named **TwoSpan_Bridge_ODB** containing the information on the nodes and elements to visualize the structure in future without using a OpenSeesPy model file.
 * a sub-folder named **Dynamic_GM1** containing the information on node displacement data to plot deformed shape.
 * a sub-folder named **ModeShapes** containing information on modeshapes and modal periods.
 * the node displacement data will be saved at closest time-steps at each 0.05 sec interval. 
   

   
