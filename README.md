# Jx_DMFT
Jx_DMFT is a software for calculating magnetic exchange parameters (Jx) from the magnetic force theorem. 
This package provides the following utilities :

* Jx calculation combined with various density functional theory (DFT) code.
* Jx calculation combined with dynamical mean-field theory (DMFT).
* Self-consistent DMFT calculation.
* Post-processing script for calculating measurable quantities (Ex. Magnon dispersion, local and momentum-dependent spectral function)

Jx_DMFT performs the above calculations by reading a non-interacting Hamiltonian through the DFTforge interface (https://github.com/KAIST-ELST/DFTforge.jl). Therefore, the list of DFT-code combinable with Jx_DMFT follows the compatibility of DFTforge.

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


    ### CTQMC of EDMFTF package installation
    Below is a very brief introduction to the installation of EDMFTF package. Therefore, it may not work in your computer environment. For a more successful installation, please visit the EDMFTF package site directly (http://hauleweb.rutgers.edu/tutorials/).
    
    Download and extract compressed EDMFTF package
    ```
    $ wget http://hauleweb.rutgers.edu/downloads/EDMFTF.tgz
    $ tar -zxvf EDMFTF.tgz
    ```
    
    Modify file configure.py to suit your computer's compiler environment
    ```
    $ vim EDMFTF-*/configure.py
    ```
    
    ,and run
    ```
    $ python setup.py
    ```
    
    If at least the following files are created in ``/bini`` folder, it is possible to run Jx_DMFT.
    ```
    ../bini/ctqmc
    ../bini/atom_d.py
    ../bini/maxent_run.py
    ```

    Finally, add the following lines to ``$ vi ~/.bashrc``
    ```
    export WIEN_DMFT_ROOT=EDMFTF-installed folder/bini
    export PYTHONPATH=$PYTHONPATH:$WIEN_DMFT_ROOT
    export SCRATCH="."
    export PATH=$WIEN_DMFT_ROOT:$PATH
    ```
    

    ### CTQMC(ComCTQMC) of COMSUITE package installation
    Below is a very brief introduction to the installation of the COMSUITE package. Therefore, it may not work in your computer environment. For a more successful installation, please visit the COMSUITE package site directly (https://github.com/comscope/comsuite).
    
    Download zip file from https://github.com/comscope/comsuite and unzip the compressed folder
    ```
    $ unzip comsuite-master.zip 
    ```
    
    Modify file arch.mk to suit your computer's compiler environment
    ```
    $ vim comsuite-master/arch.mk
    ```
    
    and, add the following two lines to ``$ vi ~/.bashrc``
    
    ```
    export COMSUITE_BIN= comsuite-installed folder/comsuite-master/bin  
    export PATH=$COMSUITE_BIN:$PATH
    ```
   
    Finally, 
    ```
    $ make all
    ```
    
    If the following files are created in ``/bin`` folder, you are ready to run Jx_DMFT.
    
    ```
    ./comsuite-master/bin/CTQMC  
    ./comsuite-master/bin/EVALSIM 
    ```
    
    
    
