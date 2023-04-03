# Power System Analysis

This Python script performs power system analysis using the Newton-Raphson method for load flow calculations.

## Overview

Power systems are complex networks consisting of generators, transmission lines, transformers, and loads. Analyzing these systems is essential for planning, operation, and maintenance of the electrical grid. Load flow analysis is a widely used method for determining the voltage, current, and power flow in an electrical network. This project implements the Newton-Raphson method, a popular iterative algorithm for solving load flow problems.

## Description

This code imports the pandas and numpy libraries and defines classes for power system components (`Barra` and `System`) and their respective functions. It reads input data from CSV files and calculates line flux, submatrices, admittance matrix, power, and updates power values. The main algorithm used is the Newton-Raphson method for load flow calculations.

### Classes

- `Barra`: Represents a bus or node in the power system. Each bus has an ID, type, voltage, angle, active and reactive power injections, and other attributes.
- `System`: Represents the entire power system and contains a list of `Barra` objects. The class contains methods for performing load flow analysis, such as calculating the admittance matrix, power mismatch, and updating voltage values.

### Input Data

The input data must be provided in two CSV files: `dadosbarra.csv` and `dadoslinha.csv`. 

`dadosbarra.csv` should contain the following columns:

- ID: Bus identifier
- Type: Bus type (slack, PV, or PQ)
- V: Initial voltage magnitude (p.u.)
- Angle: Initial voltage angle (degrees)
- P: Active power injection (p.u.)
- Q: Reactive power injection (p.u.)

`dadoslinha.csv` should contain the following columns:

- From: Sending bus ID
- To: Receiving bus ID
- R: Resistance (p.u.)
- X: Reactance (p.u.)
- B: Susceptance (p.u.)

### Output

The script will display the voltage magnitude, voltage angle, active power, and reactive power for each iteration of the Newton-Raphson method. The results can be used for various power system analyses, such as contingency analysis, voltage stability assessment, and transmission loss calculations.

## Usage

1. Prepare the input data in two CSV files: `dadosbarra.csv` and `dadoslinha.csv`.
2. Ensure that the input data is in the correct format.
3. Run the script in a Python environment.
4. The script will display the results for each iteration of the Newton-Raphson method.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change. Please make sure to update tests as appropriate.

## License

[MIT]
Copyright (c) [2023] [Matheus Domingos Ferreira da Silva]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
