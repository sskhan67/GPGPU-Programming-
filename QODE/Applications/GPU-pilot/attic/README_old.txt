0) To compile all the C files do:

    cd numerical_engines
    cmake .
    make

   If the C files are changed, 'make' needs to be re-run.
   If CMakeLists.txt are changed, refer to Yuhong's cmake email.  cmake 3.0.2 (or better) is needed.

1) Gather inputs

      1e_states.npy, 2e_states.npy, 3e_states.npy
      Z_1e.npy, Z_2e.npy, Z_3e.npy
      biorthogonal_core_matrix_dimer.npy	# Make sure to double check with Boris . . .
      biorthogonal_twobody_matrix_dimer.npy	#  . . . for semi-MO biorthogonal h and V matrices

2) Run or debug the building routines

   To run:
    python3 build_density_tensors.py 			# building density tensors (calls C code)
    python3 biorthogonal_Hamiltonian.py 		# building biorthogonal Hamiltonian (calls C code)
    python3 cmp_ham.py 					# converts Tony's Hamiltonian matrix (2e,1e,3e ... ) into Yuhong's format (2e,3e,4e,5e,6e) in order to test the code

   To debug:
    python3 test_build_density_tensors.py		# testing density tensors (calls C code)
    python3 debug_tensor.py
    python3 test_biorthogonal_Hamiltonian.py		# testing biorthogonal Hamiltonian (calls C code)
    python3 h_debug.py

3) More info:

   a. fci_index.py file is needed for indexing frozen-core FCI states.
   b. eqns/eqns.tex has derivations for building biorthogonal Hamiltonian from density tensors


