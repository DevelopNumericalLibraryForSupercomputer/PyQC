The python based quantum chemistry program.

Currently, only the RHF-based EOM-CCSD program is made and the program is tested with c++ based einsum routine.
This version moved some permuations into Hbar file, thus the einsum operations are reduced.

To run the program, c++ einsum is needed, or change es.c_einsum to numpy.einsum.

The Hbar routine and other examples will be added in future.
