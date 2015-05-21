Installation
------------

**On Linux 64-bit, you may use the provided executable**

This runs on kernel 2.6.18 and newer:
https://github.com/broadinstitute/variant-bam/src/variant

**Otherwise, you must build the executable from source**

The source code is available: https://github.com/broadinstitute/variant-bam/src

Install the dependencies:

.. code:: bash

    # Broad Institute
    use GCC-4.9
    use BamTools
    path_to_bamtools=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21	

Download and compile the code:

.. code:: bash

    #   Clone with git; easily get updates with 'git pull':
    git clone https://github.com/broadinstitute/variant-bam.git
    cd variant-bam

    cd src; 
    ./configure --with-bamtools=${path_to_bamtools}
    make               #   Compile.
    cp ../src/variant ~/bin/   #   Copy the executables wherever you like.

C++ Libraries
~~~~~~~~~~~~

To compile VariantBam, you will need a modern C++ compiler that supports
`c++0x <https://gcc.gnu.org/projects/cxx0x.html>`__ and the dependencies
listed below. I compiled successfully with gcc version 4.9.0 

`GCC, the GNU Compiler <http://gcc.gnu.org>`__

    The GNU Compiler Collection is a compiler system produced by the GNU
    Project supporting various programming languages.

