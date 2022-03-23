# vfo (Visualization For OpenSees)
Python Commands to visualize OpenSees models and components.

Install this package with,
```bash
pip install vfo
```

Upgrade the package with,
```bash
python -m pip install --upgrade vfo
```

To use this package, import the commands from **vfo** instead of *openseespy.postprocessing*. For example,

```bash
# import openseespy.postprocessing.Get_Rendering as opsplt
import vfo.vfo as vfo
```

Now, use all the openseespy visualization commands ([Here](https://openseespydoc.readthedocs.io/en/latest/src/plotcmds.html)). For example,

```bash
vfo.plot_model()
```

# USER MANUAL
A detailed list of all the latest commands will be added soon to the wiki page ([Wiki](https://github.com/u-anurag/vfo/wiki)).
