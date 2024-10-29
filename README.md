The python based quantum chemistry program.

Currently, only the RHF-based EOM-CCSD program is made and the program is tested with c++ based einsum routine.
In this version, some permuations needed in sigma vector are done when the Hbar files are generated, thus the number of einsum operation is reduced in the EOM-CCSD iteration.

This program needs c++ einsum, or one can simply change 'es.c_einsum' to 'np.einsum' for the use of the conventional einsum routine.
To run this program, type following in the command line.
python EOMCCSD.py [input file]

