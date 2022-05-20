***
# vfo (Visualization For OpenSees)
***
 
![vfo](https://github.com/u-anurag/vfo/blob/main/Doc/Shell_and_Brick3D.png)
  
  
***
[![PyPI version fury.io](https://badge.fury.io/py/vfo.svg)](https://pypi.org/project/vfo/)
[![Downloads](https://pepy.tech/badge/vfo)](https://pepy.tech/project/vfo)
[![Downloads](https://pepy.tech/badge/vfo/month)](https://pepy.tech/project/vfo)
[![Downloads](https://pepy.tech/badge/vfo/week)](https://pepy.tech/project/vfo)
***


vfo (Visualization for OpenSees) is a Python package to make your life better by helping you visualize your OpenSees models, Python or Tcl. It utilizes Matplotlib 3.0 library to plot 2D and 3D models in a dedicated interactive window. You can use click-and-hold to change the view angle and zoom the plot. The model image can be saved with the desired orientation directly from the interactive plot window. If you did not install matplotlib using Anaconda, you will have to install PyQt or PySide to enable an interactive window (Matplotlib Dependencies).

#### Animation: To save the animation movie as .mp4 file, FFmpeg codecs are required.

When using Spyder IDE and Jupyter notebook, the default setting is to produce a static, inline plot that is not interactive. To change that, write the command %matplotlib qt in the Ipython console and then execute the model plotting commands. This will produce an interactive plot in a dedicated window.

The following elements are supported:

- 2D and 3D Beam-Column Elements 
- 2D and 3D Quad Elements 
- 2D and 3D Tri Elements 
- 8 Node Brick Elements 
- Tetrahedron Elements

****
Install this package with,
```bash
pip install vfo
```

Upgrade the package with,
```bash
python -m pip install --upgrade vfo
```

To use this package, import the commands from *vfo*. For example,

```bash
import vfo.vfo as vfo
```

Now, use all the *vfo* visualization commands. For example,

```bash
vfo.plot_model()
```

# USER MANUAL

Documentation for *vfo* can be found here: ([vfo-docs](https://vfo.readthedocs.io/en/latest/index.html)). 

A detailed list of all the latest commands will be added soon to the wiki page ([Wiki](https://github.com/u-anurag/vfo/wiki)).
