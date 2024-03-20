# Tire Dynamics Repository

This repository contains scripts and data related to tire dynamics and vehicle configuration, utilizing the Pacejka Magic Formula for tire modeling.

## Structure

The repository includes the following key components:

- `tire_data/`: A directory containing tire-related data and configurations.
  - `TIR_files/`: Stores `.TIR` tire files which contain tire model parameters.
  - `config/`: Contains vehicle configuration parameters.
  - `yamls/`: Saves tire parameters in YAML format, converted from `.TIR` files using the `readTIR.py` script.
- `MagicFormula.py`: Script for computing the complete Pacejka tire model.
- `MF_lateral.py`: A lite version of the script that computes only the lateral dynamics of vehicle tires.
- `normalForce.py`: Calculates the normal force on the tire.
- `test_lateral.py`: Provides test cases for users to try out the MagicFormula and MF_lateral scripts.

## How to Use

To convert `.TIR` tire files into YAML format and compute various tire dynamics, follow these steps:

1. Place your `.TIR` tire files into the `tire_data/TIR_files/` directory.
2. Run `readTIR.py` to convert the `.TIR` files into YAML format parameters, which will be saved in the `tire_data/yamls/` directory.
3. Use `MagicFormula.py` to calculate the complete tire model dynamics, or `MF_lateral.py` for only lateral dynamics.
4. The `normalForce.py` script can be used to determine the normal force on a given tire.
5. To test and understand the implementation, run the `test_lateral.py` script.

### Prerequisites

Ensure that you have Python installed on your system as all scripts are written in Python.

## Contribution

Feel free to fork this repository and submit pull requests to contribute to the development of tire dynamics models.

## License

This project is licensed under the [Apache License](LICENSE).
