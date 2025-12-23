# BIFROST Virtual Experiments

**BIFROST Virtual Experiments** is a collection of simulation workflows and analysis notebooks for the BIFROST neutron instrumentation platform. This repository contains virtual experiment setups, instrument simulations, and analysis tools used to study instrument performance, pulse shaping, and energy resolution.

The repository is intended to support reproducible research and development related to BIFROST instrument modeling.

---

## Overview

This project includes Jupyter notebooks and supporting files for:

- Monte Carlo Particle List (MCPL) generation for BIFROST instruments
- Energy resolution studies
- Pulse shaping and tuning investigations
- Instrument configuration experiments for magnon and related simulations

The contents are primarily research-oriented and may include exploratory or legacy approaches.

---

## Repository Structure

```
.
├── BIFROST_Basis_instrument/
│   └── MCPL_generation_instrument/   # MCPL generation scripts and notebooks
├── Energy_resolution/
│   └── Old_approach/                 # Legacy energy resolution workflows
├── PSC_settings/
│   └── Magnon_instrument_sim/        # PSC and magnon simulation notebooks
├── Pulse_tuning/                     # Pulse tuning studies
├── NACALF/                           # NACALF-related scripts or data
├── Figures/                          # Generated figures (PNG, etc.)
├── .gitignore
└── README.md
```

---

## Requirements

The repository primarily relies on Python-based scientific tooling.

Recommended environment:

- Python 3.8+
- Jupyter Notebook or JupyterLab
- Common scientific Python packages (e.g. `numpy`, `scipy`, `matplotlib`)
- Instrument-specific tools such as MCPL / McStas-related libraries where applicable

Dependencies are not currently pinned in a `requirements.txt` file and may need to be installed manually depending on the notebooks used.

---

## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/KristineKrighaar/BIFROST_Virtual_Experiments.git
   cd BIFROST_Virtual_Experiments
   ```

2. Create and activate a virtual environment (optional but recommended):
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. Launch Jupyter:
   ```bash
   jupyter lab
   ```

4. Navigate to the relevant directory and open a notebook to begin.

---

## Workflows

It is often the case that the *KGS* libary is used. This is a personal libary of Kristine Krighaar and details on content can be found on. 
GitHub: [KGS](https://github.com/KristineKrighaar/KGS)

---

## Figures and Data

The `Figures/` directory contains generated plots and visualizations used by the notebooks. Raw simulation output files may not be included due to size considerations.

---

## Contributing

Contributions are welcome. If you wish to contribute:

1. Fork the repository
2. Create a feature or fix branch
3. Submit a pull request with a clear description of the changes

---

## License

No license file is currently included.  
If you intend this repository to be reused publicly, consider adding an open-source license (e.g. MIT, BSD, Apache 2.0).

---

## Author

**Kristine M. L. Krighaar**  
GitHub: [KristineKrighaar](https://github.com/KristineKrighaar)
