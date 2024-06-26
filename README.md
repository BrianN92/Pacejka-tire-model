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
2. Run `readTIR.py` from the command line, specifying the filename of the `.TIR` file you want to convert, along with optional arguments for the data and save paths if you're not using the default directories. For example:

   ```bash
   python readTIR.py -f example
   ```
   
   If your `.TIR` files are located in a different directory or you wish to save the YAML files to a different location, you can specify these paths with the `--data_path` and `--save_path` arguments, respectively:

   ```bash
   python readTIR.py example.tir --data_path '/path/to/TIR_files/' --save_path '/path/to/yamls/'
   ```

   Replace `/path/to/TIR_files/` and `/path/to/yamls/` with the actual paths where your `.TIR` files are stored and where you want to save the YAML files, respectively.

3. Use `MagicFormula.py` to calculate the complete tire model dynamics, or `MF_lateral.py` for only lateral dynamics.

4. The `normalForce.py` script can be used to determine the normal force on a given tire.

5. To test and understand the implementation, run the `test_lateral.py` script, specifying the model you wish to test. You can select between the complete tire model and the lateral tire model using the `-m` argument:

   - For the lateral tire model, use:
     ```bash
     python test_lateral.py -m lateral
     ```
   
   - For the complete tire model, use:
     ```bash
     python test_lateral.py -m MF
     ```

## References

This project draws upon the foundational works in tire and vehicle dynamics. For further reading and a deeper understanding of the concepts implemented, the following texts are highly recommended:

1. **Tire and Vehicle Dynamics** by Hans Pacejka (2005, Elsevier). This book provides an in-depth exploration of tire behavior and its impact on vehicle dynamics, offering insights into the models and theories that are central to the field.

2. **Vehicle Dynamics and Control** by Rajesh Rajamani (2011, Springer Science & Business Media). Rajamani's work delves into the principles of vehicle dynamics, control systems, and their applications in real-world scenarios, making it a valuable resource for anyone interested in the control aspect of vehicle dynamics.

These references serve as essential resources for understanding the complex interactions between tires and vehicle dynamics, and they underpin much of the theoretical framework utilized in this project.


## Contribution

Feel free to fork this repository and submit pull requests to contribute to the development of tire dynamics models.

## License

This project is licensed under the [Apache License](LICENSE).
