# PyEEDL
A tool for preparing electron transport data for use in [MCDC](https://github.com/CEMeNT-PSAAP/MCDC), derived from the Evaluated Electron Data Library ([EEDL](https://www-nds.iaea.org/epics/)).


## How to Run

### 1. Create and activate a new environment:

```bash
python3 -m venv pyeedl-env
source pyeedl-env/bin/activate
pip install endf pandas numpy h5py
git clone https://github.com/melekderman/PyEEDL.git
```

### 2. Download the EEDL data

Make sure you have the ENDF-formatted EEDL data files for individual elements (Z = 1 to 100), each stored in a separate file, and place them in the `eedl/` directory:

```bash
pip install requests bs4
cd PyEEDL/
python download_data.py
```

### 3.1 Generate the HDF5 files for raw data:

```bash
cd PyEEDL/pyeedl
python create_raw_data.py
```
### 3.2 Generate the HDF5 files for MCDC:

```bash
cd PyEEDL/pyeedl
python create_mcdc_data.py
```