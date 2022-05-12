.. include:: sub.txt

=========================
 vfo
=========================

vfo (Visualization for OpenSees) is a Python package to make your life better by helping you visualize your OpenSees models, Python or Tcl. 
It utilizes Matplotlib 3.0 library to plot 2D and 3D models in a dedicated interactive window. You can use click-and-hold to change the view angle and zoom the plot. 
The model image can be saved with the desired orientation directly from the interactive plot window. If you did not install matplotlib using Anaconda, you will have 
to install PyQt or PySide to enable an interactive window (Matplotlib Dependencies).

**Animation**: To save the animation movie as .mp4 file, `FFmpeg <https://www.ffmpeg.org/about.html>`_ codecs are required.

When using Spyder IDE and Jupyter notebook, the default setting is to produce a static, inline plot which is not
interactive. To change that, write the command **%matplotlib qt** in the Ipython console and then execute the model plotting commands. This will produce an interactive plot in a dedicated window.


Following elements are supported:

    * 2D and 3D Beam-Column Elements
    * 2D and 3D Quad Elements
    * 2D and 3D Tri Elements
    * 8 Node Brick Elements
    * Tetrahedron Elements (to be added)

The following two commands are needed to visualize the model, as shown below:

::

   #Change plot backend to Qt. ONLY if you are using an Ipython console (e.g. Spyder)
   %matplotlib qt
   
   #Change plot backend to 'Nbagg' if using in Jupyter notebook to get an interactive, inline plot.
   %matplotlib notebook
   
   # import vfo rendering module
   import vfo.vfo as vfo
   
   # render the model after defining all the nodes and elements
   vfo.plot_model()

   # plot mode shape
   vfo.plot_modeshape(3)
   

.. image:: /_static/ModelVisualization_Intro.png

Following are commands and development guide related to model visualization:

#. :doc:`createODB`
#. :doc:`saveFiberData2D`
#. :doc:`plot_model`
#. :doc:`plot_modeshape`
#. :doc:`plot_deformedshape`
#. :doc:`plot_fiberResponse2D`
#. :doc:`animate_deformedshape`
#. :doc:`animate_fiberResponse2D`
#. :doc:`plotting_OpenSeesTcl`
#. :doc:`Plotting_Development_Guide`

.. toctree::
   :maxdepth: 1
   :hidden:

   createODB
   saveFiberData2D
   plot_model
   plot_modeshape
   plot_deformedshape
   plot_fiberResponse2D
   animate_deformedshape
   animate_fiberResponse2D
   plotting_OpenSeesTcl
   Plotting_Development_Guide



