# LncLOOMv1.0 
## [Ulitsky Lab](https://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/) / [Weizmann Institute of Science](https://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/)

#### Version 1.0
#### Release Date: September 2020


Developed and maintained by [Caroline Ross](mailto:caroline-jane.ross@weizmann.ac.il) and [Igor Ulitsky](mailto:igor.ulitsky@weizmann.ac.il)

## About lncLOOM
LncLOOM is a graph-based framework that uses integer programming to identify combinations of short motifs that are 
deeply conserved in rapidly evolving sequences. This version is implemented in Python 2 and is supported on Linux/Unix-based
systems. 

## Getting Started (Install with pip)
*If you're unable to install LncLOOM with pip, please see examples in the Troubleshooting section at the bottom of this file on how to run LncLOOM from within  the LncLOOM_v1 directory.

1. Download the lncLOOM repository.

   `git clone https://github.com/lncLOOM/LncLOOMv1.git`
   
   `cd LncLOOMv1`

2. Install Python 2 (If needed)
   LncLOOM is supported on Linux/Unix-based systems. It is run via the command line. 
   Python 2 is available for download at:
   <A HREF="https://www.python.org/downloads/">here</A>

3. Install LncLOOM as an executable using pip.

   * Firstly ensure that [pip](https://pypi.org/project/pip/) is installed:
   
     `sudo apt-get install python-pip`

      if you are using macOS:
      
     `sudo easy_install pip`
     
     
    * Install LncLOOMv1 using pip (the following command ensures that it is setup to run with python2) 
    
       `python2 -m pip install --user -e ./LncLOOM_v1`

     
4. Add LncLOOM to your $PATH. 
   
    pip creates a LncLOOM executable. Depending on your OS, this executable will be saved to certain directory, which needs to be added to your $PATH:
      
    * For Linux systems, LncLOOM will be saved in ~/.local/bin/ (or /home/<you>/.local/bin)
      
         `export PATH="~/.local/bin:$PATH"`    
          
      
    * For macOS, LncLOOM will be saved in /Users/<you>/Library/Python/<version>/bin
      eg: /Users/Mac/Library/Python/2.7/bin
      
         `export PATH="/Users/Mac/Library/Python/2.7/bin:$PATH"`
   
    *Note: the above paths to LncLOOM may vary depending on your directories
    
    
5. Install LncLOOM dependencies

   LncLOOM requires several packages to be installed. Most of these would have already been installed when you installed LncLOOMv1 
   (see last section of this page for a list of these packages)
   However the following additional programs must be installed individually:
   
   * BLASTN 
   
   `sudo apt-get install ncbi-blast+`
   
   * [Mafft](https://mafft.cbrc.jp/alignment/software/)  
     To download and setup follow the steps given [here](https://mafft.cbrc.jp/alignment/software/linux.html)
     


6. Set paths to genome files and eCLIP data that LncLOOM will use for annotations and generation of a custom track for the UCSC Genome Browser 
   * In the LncLOOM_v1/LncLOOM_v1/src/ directory there is a file called `for_eclip_annotation.txt`. This file tells LncLOOM where to find data needed for annotations. The file looks as follows:
     ```
           Query Layer: 1
           Blat: src/hg19.fa
           eCLIP: Data 1: src/eCLIP_narrowPeakApr2019/
     ```
     Currently the paths have been set to use data that is located in the LncLOOM_v1/LncLOOM_v1/src/ folder. However, these files are too large to be stored on GitHub and need to be downloaded from [hg19.fa](https://drive.google.com/file/d/1ZyXCX1o7S4g0Ad_6wZU9esRNVzeqiNCc/view?usp=sharing) and [eCLIP_narrowPeakApr2019](https://drive.google.com/file/d/1PaU5kJvfC26fENF6E0H-fwLTimlW-e-8/view?usp=sharing)
     
     To use this data, download and extract the files into the LncLOOM_v1/LncLOOM_v1/src/ folder.
     The eCLIP data consists of BigBed files retrieved from [ENCODE](https://www.encodeproject.org/eclip/) in 2019.
          
     `tar xvzf eCLIP_narrowPeakApr2019.tar.gz`
     
     `mv eCLIP_narrowPeakApr2019 LncLOOMv2/LncLOOMv2/src/`
     
     `tar xvzf hg19.tgz`
     
     `mv hg19.fa LncLOOMv2/LncLOOMv2/src/`
     
     Alternatively, if you have your own data you can update these paths in `for_eclip_annotation.txt` to the full paths to your genome file and eCLIP data. For example:
     ```
           Query Layer: 1
           Blat: /home/MySpace/MyGenomeFiles/hg19.fa
           eCLIP: Data 1: /home/MySpace/My_eCLIP_Data/
     ```
     
     To annotate motifs found with eCLIP data specified in `for_eclip_annotation.txt` use the `--eclip` option when running LncLOOM.
     
     - Explanation:
       * The query layer specifies which sequence you would like annotate. By default this will be the top sequence (layer 1) in your input file.
              Note that LncLOOM always sets the first sequence in your file to the top sequence, but may reorder the other sequences to improve motif
              discovery. To retain your original order of sequences use the `--inputorder` command is used. 

       * Blat: specifies the full path to a genome file

       * eCLIP: specifies the full paths to eCLIP data. Note: you can add multiple paths as follows:

         ```
         Query Layer: 1
         Blat: <specify path to genome fasta file>
         eCLIP: Data 1: <specify path to eCLIP data>
         eCLIP: Data 2: <specify path to eCLIP data>
         eCLIP: Data 3: <specify path to eCLIP data>
         ```

         Alternatively you can upload a bedfile instead of running Blat
         ```
         Query Layer: 1
         Bed: <specify path to bedfile>
         eCLIP: Data 1: <specify path to eCLIP data>
         ```

In the LncLOOM_v1/LncLOOM_v1/src/ directory there is also a file: `for_track_output.txt`. Similar to the `for_eclip_annotation.txt`, this file tells LncLOOM where to find a genome file so that a custom track of conserved motifs can be generated. The paths have been initiated to find `hg19.fa` in LncLOOM_v1/LncLOOM_v1/src. 
Note that you can specify a different layer and genome to what is specified in  `for_eclip_annotation.txt`.

         To generate a custom track use the `--track` option.
         ```
         Query Layer: 1 
         Blat: <specify path to genome fasta file>
         
         
7. Make sure that the blat executable has the correct executable permissions:

   `chmod 755 LncLOOM_v1/LncLOOM_v1/src/blat`
         
8. OPTIONAL: Install the [Gurobi Solver](https://www.gurobi.com/) - although not required it allows much faster performance on larger datasets
   There are two possible ways to install Gurobi:
   - Option 1: Install through [Anaconda](https://www.gurobi.com/gurobi-and-anaconda-for-linux/).  
     - If needed download and install [Anaconda](https://www.anaconda.com/products/individual)  
     - Add the gurobi channel to the Ananconda search list   
       `conda config --add channels http://conda.anaconda.org/gurobi`
     - install gurobi  
       `conda install gurobi`
     - Initialise Gurobi License 
       - A free academic license can be obtained from: [https://www.gurobi.com/downloads/end-user-license-agreement-academic/]
       - First register an account
       - Verify your account from a link sent to your email, this will take you to a home page
       - Click on Licenses (you may be askd to login again)
       - On the top Navigation bar, select Academia... and Licenses
       - Click on a link: Free Academic License page,this will issue you a license
       - Scroll to the bottom of the page to Installation: you will see a command similar to this, but specific to your key: `grbgetkey YOUR_KEY`
       - Copy and run this command in your terminal  
   - Option 2:  
     - Download [Gurobi](https://www.gurobi.com/downloads/gurobi-software/)

     - Once you have downloaded your version of Gurobi copy the folder to `/opt`

       `sudo cp -r gurobi9.0.2_linux64.tar.gz /opt`
     - Extract the file into `/opt`

       ```
       cd /opt/
       sudo tar xvfz gurobi9.0.2_linux64.tar.gz
       ```

      - Set environment variables  

        Users of the bash shell should add the following lines to their .bashrc files:
        ```
        export GUROBI_HOME="/opt/gurobi902/linux64"
        export PATH="${PATH}:${GUROBI_HOME}/bin"
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
        ```
        
        Users of the csh shell should add the following lines to their .cshrc files:
        ```
        setenv GUROBI_HOME /opt/gurobi902/linux64
        setenv PATH ${PATH}:${GUROBI_HOME}/bin
        setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib
        ```
        
        If `LD_LIBRARY_PATH` is not already set, use the following instead:
        
        `export LD_LIBRARY_PATH="${GUROBI_HOME}/lib"`
        or
        `setenv LD_LIBRARY_PATH ${GUROBI_HOME}/lib`


     - Initialise Gurobi License 
       - A free academic license can be obtained from: https://www.gurobi.com/downloads/end-user-license-agreement-academic/
       - First register an account
       - Verify your account from a link sent to your email, this will take you to a home page
       - Click on Licenses (you may be askd to login again)
       - On the top Navigation bar, select Academia... and Licenses
       - Click on a link: Free Academic License page,this will issue you a license
       - Scroll to the bottom of the page to Installation: you will a command similar to this, but specific to your key:
            `grbgetkey YOUR_KEY`
       - copy and run this command in your terminal
         
## Running LncLOOM

* Basic command to run LncLOOM which invokes all default options:
  ```
   LncLOOM --fasta <path to file of sequences>
  ```
* To save output in a specified directory use the `--pname` command:
  ```
   LncLOOM --fasta <path to file of sequences> --targetscan --pname <name of directory>
  ```
* To annotate motifs with conserved miRNA binding sites from [TargetScan](http://www.targetscan.org/), invoke the `--targetscan` option:
  ```
   LncLOOM --fasta <path to file of sequences> --pname <name of directory> --targetscan
  ```

* To annotate motifs based on eCLIP data according to the parameters set in `src/for_eclip_annotation.txt`, use the `--eclip` option:
  ```
   LncLOOM --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip
  ```
* To perform an empirical statistical analysis, run multiple iterations on random datasets generated by LncLOOM:
  ```
   LncLOOM --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip --iterations 100
  ```
* To speed things up you can run the statistical iterations in parallel
  ```
  LncLOOM --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip --iterations 100 --multiprocess 10
  ```
  
* To generate a custom track of your conserved motifs,coloured by conservation, that can be viewed in Genome Browser, use the `--track` command:
  Output will be in bedfile format
  ```
  LncLOOM --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip --iterations 100 --multiprocess 6 --track
  ```
## Defintions of Statistical values

Definitions and troubleshooting tips for calculating motif significance are also given in [Definitions.html](https://github.com/lncLOOM/lncLOOM/edit/master/Definitions.html):

P(i): Probability of finding the exact motif, at the same depth, in a random set of sequences that have the same percentage identities as the input sequences

E(i): Probability of finding any combination of the same number of motifs of the same length, or longer, at the same depth, in a random set of sequences that have the same percentage identities as the input sequences

P(r): Probability of finding the exact motif, at the same depth, in a random set of sequences that have the dinucleotide composition as the input sequences

E(r): Probability of finding any combination of the same number of motifs of the same length, or longer, at the same depth, in a random set of sequences that have the same dinucleotide composition as the input sequences


## All LncLOOM Options
LncLOOM has several options:

### Required arguments
  
* `--fasta` path to fasta file



### Optional arguments 

* `--similarity` type is a float, removes all sequences that have identity higher than the specification. Default=100.00

* `--mindepth` LncLOOM will only find motifs that are conserved to this depth. Default=2

* `--maxedges` Maximum number of edges allowed in any graph, default = 1200, if number edges exceeds this then no solution is returned.

* `--solver` *gurobi* or *cbc* (default = cbc) - gurobi is much faster.

* `--startw` LncLOOM will begin scanning for motifs of this length and then iteratively decrease to a length = stopw. Default=15

* `--stopw` Shortest length of kmer that LncLOOM will search for. Default=6

* `--outdir` Directory for output. Default: current working directory.

* `--pname` Name of project, folder is created within outdir and all files are written to this folder. Default: LincMotif

* `--select` Name of sequence header: LncLOOM will print out all motifs found to the depth of this sequence. Default is to print out deepest kmers found 

* `--prune` LncLOOM will emove all kmers that appear more than this number of times in a single sequence. Default = 8

* `--iterations` Number of iterations to run for compting P values. Default = 0

* `--multiprocess`: LncLOOM will divide the iterations into this number of processes. Default = 4

* `--tol5`: Tolerance step from median postion of first and last nodes to determine exclusion from 5' extension graphs. Default = 0.50
    
* `--tol3`: Tolerance step from median postion of first and last nodes to determine exclusion from 3' extension graphs. Default = 0.50


### BOOLEAN OPTIONS

The following are boolean options (all defaults are false, by simple typing --option, it will be set to true)

* `--targetscan` boolean (default=False). If set to true motifs will be annotated with conserved miRNA sites, from data obtatined from TargetScan.
* `--eclip`: boolean (default=False). If set to true, motifs in the specified layer (see `src/for_eclip_annotation.txt`) will be annotated using eCLIP data.
* `--inputorder` if set to true, the sequences are not reordered based on homology and the graph is built according to fasta file order of sequences
* `--track` generates a custom track of conserved motifs, coloured by conservation, that can be viewed in Genome Browser (https://genome.ucsc.edu/)
* `--hspblast` boolean (default=False), if hspblast is true then the High Scoring Segment Pairs are used as constraints to build graph and edges
* `--noconstraints` shallower simple paths will NOT be used as constraints - much slower when used on larger datasets
* `--newcolours` generates output results with a new random colour scheme
* `--shorttandem` Excludes tandem repeats (complex paths) found in iterations of larger k


## Command Line Examples

More examples of commands:

   ```
   LncLOOM --fasta Chaserr.fas --name Chaserr  --startw 10  --solver gurobi --iterations 100 --multiprocess 10 --eclip --targetscan --track
   ```
   ```
   LncLOOM --fasta Chaserr.fas --name Chaserr  --startw 10  --solver gurobi --iterations 100 --eclip --targetscan --track --tol5 0.1 --tol3 0.5
   ```
   ```
   LncLOOM --fasta Chaserr.fas --pname Chaserr  --startw 10  --solver gurobi --eclip --targetscan --inputorder`
   ```
   ```
   LncLOOM --fasta Chaserr.fas --pname Chaserr  --startw 10  --solver gurobi --eclip --targetscan --inputorder --newcolours`
   ```


## Troubleshooting (If installation with pip was not successful)
   * You may need to [upgrade](https://pip.pypa.io/en/stable/installing/#upgrading-pip) to a newer version of pip
      
   Note: The following packages should have been automatically installed in step 3. However, if this was unsuccessful each package must be installed individually:
   * [networkx](https://networkx.github.io/)  
    `pip install networkx`
    
   * [PulP](https://github.com/pulp/pulp)  
    `pip install pulp`
     You may need to run the [pulp tests](https://pypi.org/project/PuLP/) to ensure that the Linear Programming solvers are available to pulp:
     
     `sudo pulptest`
     
     If PULP_CBC_CMD is unavailable, reinstall pulp using the following command, and then rerun the pulptest:
     
     `python2 -m pip install -U git+https://github.com/coin-or/pulp`
     
     If you are using Gurobi, ensure that you have [configured](https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html) the solver correctly with PulP
   
   * [NumPy](https://numpy.org/)  
    `pip install numpy`
    
   * [pyBigWig](https://github.com/deeptools/pyBigWig)  
    `pip install pyBigWig`
    
   * [Biopython](https://biopython.org/) (Version 1.76 is the last release to support python 2.7, later versions are supported by Python 3)  
    `pip install biopython==1.76`

* Run lncLOOM from within the downloaded LncLOOM_v1 directory. This is important as LncLOOM uses files stored in the src folder.
* Change your directory to LncLOOM_v1/LncLOOM_v1:
  ```
   cd ./lncLOOM/LncLOOM_v1/LncLOOM_v1
  ```
* Basic command to run LncLOOM which invokes all default options:
  ```
   python LncLOOM.py --fasta <path to file of sequences>
  ```
* To save output in a specified directory use the `--pname` command:
  ```
   python LncLOOM.py --fasta <path to file of sequences> --targetscan --pname <name of directory>
  ```
* To annotate motifs with conserved miRNA binding sites from [TargetScan](http://www.targetscan.org/), invoke the `--targetscan` option:
  ```
   python LncLOOM.py --fasta <path to file of sequences> --pname <name of directory> --targetscan
  ```

* To annotate motifs based on eCLIP data according to the parameters set in `src/for_eclip_annotation.txt`, use the `--eclip` option:
  ```
   python LncLOOM.py --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip
  ```
* To perform an empirical statistical analysis, run multiple iterations on random datasets generated by LncLOOM:
  ```
   python LncLOOM.py --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip --iterations 100
  ```
* To speed things up you can run the statistical iterations in parallel
  ```
  python LncLOOM.py --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip --iterations 100 --multiprocess 10
  ```
  
* To generate a custom track of your conserved motifs,coloured by conservation, that can be viewed in Genome Browser, use the `--track` command:
  Output will be in bedfile format
  ```
  python LncLOOM.py --fasta <path to file of sequences> --pname <name of directory> --targetscan --eclip --iterations 100 --multiprocess 6 --track
  ```

