# Power System Analysis

This project provides a set of tools for power system analysis using the Newton-Raphson and Gauss-Seidel methods to solve the power flow problem. It relies on two input databases, one for bus data and another for line data, provided as CSV files.

## Databases

The project uses two databases in the form of CSV files to store the input data for the power system analysis:

1. `dadosbarra.csv`: Contains the bus data, including the bus ID, voltage magnitude, voltage angle, active power, reactive power, and the type of bus (generator, load, or reference).

2. `dadoslinha.csv`: Contains the line data, including the IDs of the two connected buses, line resistance, line reactance, and line charging susceptance.

Both files should be placed in the same directory as the `power_system_analysis.py` file.

## Getting Started

These instructions will help you set up the project and run it on your local machine for development and testing purposes.

### Prerequisites

To use this project, you need to have the following software installed:

- Python 3.7 or later
- NumPy
- Pandas

### Installing

To install the required libraries, run the following command:

\```sh
pip install numpy pandas
\```

## Usage

The main code can be found in `power_system_analysis.py`. The project is structured with the following classes and functions:

1. Class `Barra`: Represents a bus object in the power system. Each object contains information about the bus ID, voltage magnitude, voltage angle, active power, reactive power, and the type of bus.

2. Class `System`: Represents the entire power system. This class contains a list of buses, a list of lines, and the admittance matrix of the system.

3. Function `read_csv_data`: Reads the input data from two CSV files, `dadosbarra.csv` for bus data, and `dadoslinha.csv` for line data. It returns two Pandas DataFrames containing the data.

4. Function `create_system`: Takes the bus and line DataFrames and creates a `System` object containing a list of `Barra` objects and a list of lines.

5. Function `calculate_admittance_matrix`: Takes the system object and calculates the admittance matrix based on the bus and line data.

6. Function `calculate_power`: Calculates the active and reactive power of each bus in the system.

7. Function `newton_raphson_method`: Implements the Newton-Raphson method to solve the power flow problem. It iteratively updates the voltage magnitudes and angles until the system converges.

8. Function `display_results`: Displays the power flow results, including voltage magnitudes, voltage angles, active power, and reactive power.

9. Function `update_power_values`: Updates the active and reactive power values of the buses in the system after each iteration of the Newton-Raphson method.

10. Function `gauss_seidel_initial_values`: Initializes the voltage magnitudes and angles for the Gauss-Seidel method.

### Running the program

To run the program, execute the following command:

\```sh
python power_system_analysis.py
\```

The program will read the input data from `dadosbarra.csv` and `dadoslinha.csv` and display the results in the terminal.

## Contributing

We welcome contributions to improve the project. Please follow these steps to contribute:

1. Fork the repository.
2. Create a new branch with a descriptive name.
3. Make your changes in the new branch.
4. Submit a pull request to the main branch.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.


