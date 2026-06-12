# Stress-constrained Variable Thickness Topology Optimization with deal.II

This repository contains a modular C++ framework for 2D MPI-parallel Topology Optimization, built upon the [deal.II](https://www.dealii.org/) finite element library.

The code supports various methods and physics, including:

* **Methods:** SIMP, Multi-thickness, and Variable Thickness implementations.
* **Physics:** Linear and Nonlinear elasticity.
* **Regularization:** PDE-based and Heaviside filters.
* **Mesh:** Adaptive Mesh Refinement and Coarsening.

## 1. Installation

### Prerequisites

To compile and run this code, you must have docker installed. The code uses deal.II 9.4.0 FEM library, hence the docker image of deal.II must be ran.

### Cloning the Repository
Clone the code from the GitHub repository to your local machine:

```bash
git clone https://github.com/gabstank/Stress-constrained-VTTO.git
cd Stress-constrained-VTTO
```

### Docker image

This code is tested with deal.II library version 9.4.0. Therefore, you need to pull the following docker image:

```bash
docker pull dealii/dealii:v9.4.0-jammy
```

## 2. Compilation

Make sure the docker image is active and you are within the docker workspace. We use CMake to configure the build. You should create a separate build directory to keep the source tree clean.

1.  **Create a build directory:**

    ```bash
    mkdir build
    cd build
    ```

2.  **Configure the project:**
    Choose your build type based on your needs.

    * **DebugRelease:** High optimization levels but keeps debug symbols/assertions. Best for running simulations while maintaining safety checks.
    
        ```bash
        cmake -DCMAKE_BUILD_TYPE=DebugRelease ..
        ```

    * **Debug:** Unoptimized, full assertions. Use this if you are developing new features or fixing crashes.
    
        ```bash
        cmake -DCMAKE_BUILD_TYPE=Debug ..
        ```

    * **Release (Recommended):** Maximum optimization, no safety checks.
        ```bash
        cmake -DCMAKE_BUILD_TYPE=Release ..
        ```

3.  **Compile the code:**
    Use `make` with the `-j` flag to compile in parallel (replace `4` with your number of cores):
    
    ```bash
    make -j4
    ```

## 3. Running the Code

The executable `topopt` is generated in the `build/` directory.

### Directory Structure for Runs
It is recommended to run the code from a subdirectory (e.g., `examples/`) where your `.prm` (parameter) files are located, rather than from the build folder.

### Example Command
Assuming you are in the root directory `topopt/`:

1.  Navigate to the project folder:

    ```bash
    cd project
    ```

2.  Run the code using MPI (e.g., on 8 processors) pointing to the build directory:

    ```bash
    mpirun -np 8 ../build/topopt test.prm
    ```

### Output
* **Console:** Iteration progress, objective function values, and solver convergence details.
* **Files:** Main visualization file `<name>-output.pvd` will located in the `results/` folder. If the `results/` is not created manually, the optimization might cause an error. The file `<name>-convergence.csv`, generated after the optimization has finished, contains the important iteration data like the objective, design change etc. The individual iteration files `.vtu` are in the `results/VTUs` - you do not need to worry about them.

## 4. Repository Structure

* `src/`: Contains all `.h` header and `.cc` implementation files.
	* `TopOpt.h/.cc`: Main topology optimization class handling the complete optimization loop
	* `BVP/`: Boundary value problems files, for instance, linear elasticity FEM.
	* `DesignField/`: Design field related files - contain the density fields and penalization techniques.
	* `Filters/`: Density filtering implementations.
	* `Materials/`: Material law classes.
	* `Mesh/`: Mesh related classes - generation of the geometry, adaptive meshing.
	* `Optimizer/`: Optimizer classes.
	* `Responses/`: Computation of the responses and their sensitivities, e.g. compliance, volume.
	* `Utilities/`: Various utility classes: parameter handler, macros, conversion functions, global objects and more.
* `examples/`: Example parameter files.

## Publications

Relevant publications

1. *Stankiewicz, G., Dev, C., & Steinmann, P. (2026). Deblurring structural edges in variable thickness topology optimization via density-gradient-informed projection. Structural and Multidisciplinary Optimization.*
2. *Stankiewicz, G., Dev, C., & Steinmann, P. (2026). Deblurring structural edges in variable thickness topology optimization via density-gradient-informed projection. Structural and Multidisciplinary Optimization.*
