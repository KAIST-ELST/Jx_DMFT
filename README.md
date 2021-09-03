# Jx_DMFT
Jx_DMFT is a software for calculating magnetic exchange parameters (Jx) from the magnetic force theorem. 
This package provides the following utilities :

* Jx calculation combined with various density functional theory (DFT) code.
* Jx calculation combined with dynamical mean-field theory (DMFT).
* Self-consistent DMFT calculation.

Jx_DMFT performs the above calculations by reading a non-interacting Hamiltonian through the DFTforge interface. Therefore, the list of DFT-code combinable with Jx_DMFT follows the compatibility of DFTforge (https://github.com/KAIST-ELST/DFTforge.jl).

# Prerequisite
* DFTforge pakage of Julia language
* DMFT solver (now two CTQMC solver are compatible : CTQMC of COMSUITE package (https://github.com/comscope/comsuite), and CTQMC of EDMFTF (http://hauleweb.rutgers.edu/tutorials/) are available.)

  ## Installation of prerequisite
    ### DFTforge installation
    Download the latest compressed Julia Binaries from (https://julialang.org/downloads/) and extract the file. For Linux system,
    ```
    $ wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.2-linux-x86_64.tar.gz
    $ tar -xvzf julia-1.6.2-linux-x86_64.tar.gz
    ```
    and, copy the extracted folder to /opt and make a symbolic link
    ```
    $ sudo cp -r Julia-1.6.2  /opt/
    $ sudo ln -s /opt/Julia-1.6.2/bin/Julia  /usr/local/bin/Julia
    ```
    Then, you can start Julia by re-opening the terminal and typing
    ```
    $ julia
    ```
    Finally, you can install DFTforge by using Pkg.add command in Julia
    ```
    $ julia
       _       _ _(_)_     |  Documentation: https://docs.julialang.org
      (_)     | (_) (_)    |
       _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 1.6.2 (2021-07-14)
     _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
    |__/                   |
    julia> import Pkg; Pkg.add("DFTforge")
    ```
    
    ### CTQMC(ComCTQMC) of COMSUITE package installation
    Below is a very brief introduction to the installation of the COMSUITE package. Therefore, it may not work in your computer environment. For a more successful installation, please visit the COMSUITE package site directly (https://github.com/comscope/comsuite).
    ### CTQMC of EDMFTF package installation
    Below is a very brief introduction to the installation of EDMFTF package. Therefore, it may not work in your computer environment. For a more successful installation, please visit the COMSUITE package site directly (http://hauleweb.rutgers.edu/tutorials/).
