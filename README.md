# GravityToolsLight

Welcome to `GravityToolsLight`, an open-source Julia library designed to facilitate the precise timing analysis of pulsars with a focus on gravitational studies. This library serves as a comprehensive toolkit for astronomers and physicists, providing an extensible framework for running and managing the `tempo` and `tempo2` programs, exploring Post-Keplerian parameters, and performing advanced grid-based analyses.

## Features

- **GeneralTempoFramework**: A robust wrapper for running the `tempo` and `tempo2` programs, offering functionalities to manage parameter files, execute timing analysis, and iterate over model fitting procedures.
- **GeneralPKFramework**: A set of tools for working with Post-Keplerian parameters to conduct theoretical tests of gravity without direct TOA data, comparing measured values against specific gravitational theories such as General Relativity and alternative frameworks.
- **AdaptiveRefinement2DGrid**: Algorithms capable of performing adaptive refinement on a two-dimensional grid of parameters, enabling precise exploration of parameter spaces both in single-threaded and parallel computing environments.
- **Combined Frameworks**: Integration of the above frameworks to facilitate complex tasks such as grid-based `tempo` execution, PK parameter grid searches, and enhanced convergence methods combining both Tempo and PK approaches.

## Getting Started

Please refer to the [Documentation](#) for a detailed guide on installing and using `GravityToolsLight`. Examples and tutorials can also be found in the `/examples` directory to help you get started with common tasks.

## Installation

`GravityToolsLight` can be installed using the Julia package manager. From the Julia REPL, type the following:

```julia
using Pkg
Pkg.add("GravityToolsLight")
```

## Contributing

Contributions are more than welcome! If you have suggestions for improvements or want to contribute code, please feel free to fork the repository, make your changes, and submit a pull request.

## License

GravityToolsLight is licensed under the MIT License. Feel free to use it, contribute, and spread the word!

## Acknowledgments

GravityToolsLight is being actively developed and maintained by [Your Name] and contributors from the scientific community. We thank all those who have offered valuable insights, contributed code, and used the library in their research.


