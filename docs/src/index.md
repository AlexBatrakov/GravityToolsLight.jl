# GravityToolsLight Documentation

```@contents
"# GravityToolsLight Documentation

Welcome to the documentation for `GravityToolsLight`, a Julia library aimed at facilitating the analysis of pulsar timing for gravitational studies. This library provides a set of tools for running `tempo` and `tempo2`, handling Post-Keplerian parameters, and performing grid-based analyses with adaptive refinement.

## Introduction

`GravityToolsLight` is designed to be a comprehensive toolkit for astronomers and physicists who need to perform precise timing analyses of pulsars. It offers a range of functionalities from basic timing to advanced model fitting and parameter estimation.

## Installation

To install `GravityToolsLight`, use the Julia package manager:

```
using Pkg
Pkg.add("GravityToolsLight")
```


## Getting Started

Here's how to get started with `GravityToolsLight`:

using GravityToolsLight

## Example of basic usage


## Features

### GeneralTempoFramework

- Manage parameter files for `tempo` and `tempo2`.
- Execute timing analysis with comprehensive fitting procedures.

### GeneralPKFramework

- Work with Post-Keplerian parameters for theoretical gravity tests.
- Compare measured values with predictions from various gravitational theories.

### AdaptiveRefinement2DGrid

- Perform adaptive refinement on a two-dimensional grid of parameters.
- Utilize single-threaded or parallel computing environments for grid exploration.

### Combined Frameworks

- Integrate Tempo and PK methods for enhanced analysis workflows.
- Conduct grid searches and optimize parameter estimation.

## Usage Examples

Here are some common usage scenarios for `GravityToolsLight`:

### Running Tempo Analysis

### Code example for running a Tempo analysis


### Post-Keplerian Parameter Estimation

### Code example for Post-Keplerian parameter estimation

### Code example for adaptive grid refinement


## Contributing

We welcome contributions to `GravityToolsLight`! If you'd like to contribute, please fork the repository, make your changes, and submit a pull request.

## License

`GravityToolsLight` is released under the MIT License. See the [LICENSE](LICENSE) file for more details."

