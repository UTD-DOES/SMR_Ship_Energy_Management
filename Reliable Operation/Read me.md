# Reliable Operation of Dynamic Model

This repository provides all necessary files and scripts for the dynamic modeling and automatic simulation of a system integrating a Battery Energy Storage System (BESS) and SMR with dynamic propulsion motor load behavior.

## Repository Content

- **Grid_v33_B0.sav**  
  Network configuration file for the system.

- **Modeling_v2.dyr**  
  Dynamic models including the battery and SMR.

- **Operation_B2025.out**  
  Simulation output results file.

- **Operation_Script.py**  
  Python script that **automatically loads** all necessary files (`.sav`, `.dyr`, and `P_motor6.csv`) and runs the simulation.

- **P_motor6.csv**  
  CSV file providing the dynamic motor load profile.

## Overview

The provided script (`Operation_Script.py`) handles:
- Loading the **network model (.sav)**.
- Loading the **dynamic model (.dyr)**.
- Importing the **motor load profile (P_motor6.csv)**.
- Running the simulation automatically without manual file setup.

## Requirements

- Python 3.x
- Power system simulation software (e.g., **PSSE**) supporting `.sav` and `.dyr` formats
- Python libraries:
  - `pandas`
  - `numpy`
  - Simulation-specific APIs (e.g., `psspy` for PSSE)

## How to Run

1. Clone/download this repository.
2. Open `Operation_Script.py`.
3. Ensure paths to `.sav`, `.dyr`, and `.csv` are correct (if running outside the same directory).
4. Run `Operation_Script.py`.  
   It will **automatically load** the system files and start the simulation.
5. Review the simulation results in `Operation_B2025.out`.

## Applications

- Dynamic stability analysis
- Battery Grid integration studies
- Motor load impact assessments
- Reliable operation of marine energy systems under dynamic conditions

