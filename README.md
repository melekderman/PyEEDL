# PyEEDL
A tool for preparing electron transport data for use in [MCDC](https://github.com/CEMeNT-PSAAP/MCDC), derived from the Evaluated Electron Data Library ([EEDL](https://www-nds.iaea.org/epics/)).


## How to Run: 

Before running, make sure you have the ENDF-formatted EEDL data files for individual elements (Z = 1 to 100), each stored in a separate file, and place them in the `/eedl` directory.

```bash
python3 -m venv pyeedl-env
source pyeedl-env/bin/activate
pip install endf
git clone https://github.com/melekderman/PyEEDL.git
cd PyEEDL/pyeedl
python run.py
```
