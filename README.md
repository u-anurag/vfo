# openseespyvis
A developmenet branch of OpenSeesPy Visualization commands. This can be used while waiting for the next release of OpenSeesPy.
It is advisable to use the latest release of OpenSeesPy (once released) for stable functions.

Install this package with,
```bash
pip install openseespyvis
```

Upgrade the package with,
```bash
python -m pip install --upgrade openseespyvis
```

To use this package, import the commands from **openseespyvis** instead of *openseespy.postprocessing*. For example,

```bash
# import openseespy.postprocessing.Get_Rendering as opsplt
import openseespyvis.Get_Rendering as opsplt
```

Now, use all the openseespy visualization commands ([Here](https://openseespydoc.readthedocs.io/en/latest/src/plotcmds.html)). For example,

```bash
opsplt.plot_model()
```

# USER MANUAL
Check out the ([Wiki](https://github.com/u-anurag/openseespyvis/wiki)) for all the commands available in openseespyvis.
