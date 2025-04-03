# MiniSPICE
[简体中文](https://github.com/rid-sun/MiniSPICE/readme.md)/English

This project evolved from a parser and was rewritten using C++ features (inheritance, encapsulation, abstraction, and polymorphism) to restructure the overall architecture. It provides an excellent demo for <u>**easily understanding the SPICE development process**</u>. 😊😊😊

> Main references:
> * https://www.sfu.ca/~ljilja/cnl/projects/Homotopy/  
> * https://ngspice.sourceforge.io/

[![Watch the video](pic/demo.png)](https://www.bilibili.com/video/BV1thP5eXEEn/)

---------------------------------------------------------------
## Key features
- The overall architecture follows a **batch processing model**, where functionalities are specified through netlist statements.
- Constructs KCL equations **node-centric**, utilizing a **Jacobian matrix filling approach** instead of a DC companion model.
- Supports **SUPERLU_DIST** and **KLU** as **optional** solvers, enabling both distributed parallel and serial solving.
- Implements **various** matrix formats with **only one conversion** throughout the entire process, minimizing the overhead of introducing different solvers.
- Supports **multiple DC analysis algorithms**, including the traditional NR iteration and pseudo-transient analysis.
- Provides **graphical visualization** of analysis results.

-------------------------------------------------------------

## TODO
### Overall Architecture
- [ ] Interactive design
### Functional Modules
- [ ] Support for homotopy analysis using the arc-length method
- [ ] Support for transient analysis in the time domain
### Solver Module
- [ ] Development of **MiniSPICE-LU** with CUDA-accelerated left-looking LU decomposition

-----------------------------------------------------

## Build and Run Instructions
### Environment Dependencies
* Intel HPC Toolkit
* Cmake
* Python
* METIS/ParMETIS

### Build
To ensure environment consistency, it is recommended to build the project using Docker. The project has been fully tested on Windows, requiring approximately **15GB of image storage space and 410s for project compilation**.
> Development environment: `Windows + Docker + VSCode`

On Linux, the project can also be built using `OpenMPI + GCC`, but modifications to the dependencies in [`src/solver/CMakeLists.txt`](./src/solver/CMakeLists.txt) are required to ensure compatibility.

`Windows + Docker` Build Process
1. Install `Docker`: Download from the official website and configure the image storage path as needed.
2. Clone the project repository:
    ```bash
    git clone https://github.com/rid-sun/MiniSPICE.git
    ```
3. Start `Docker` and execute the following commands in order:
    ```bash
    cd MiniSPICE
    docker build . --file Dockerfile --tag minispice
    ------waiting------
    docker run --shm-size=4gb --mount type=bind,source="%CD%",target=/root/minispice -it minispice
    cd /root/minispice/3rd_lib
    bash install.sh
    ------waiting------
    cd ../
    bash build.sh
    ```
4. Build complete

### Run
```bash
cd bin
mpirun -np num ./minispice -f ../testcase/testcase1/Netlist1.txt -o netlist1
```
where `num` is the number of processes. By default, the project uses `KLU` as the solver. To enable `SuperLU_DIST`, add `-DUSE_SUPERLU=ON` during the build process.
