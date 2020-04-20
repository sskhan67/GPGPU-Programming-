#    (C) Copyright 2018 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#
import sys
import copy
import signal
import numpy as np
import traceback
from time import time
from qode.coupled_oscillators import oscillator_system as osys
from qode.coupled_oscillators import hamiltonian_class
from qode.coupled_oscillators.analytical import analytical_solution
import simple_fast
from testLanczos import experiment, build_vector, build_matrix



def get_mol_obj(
                num_oscillator,
                num_molecule,
                positions
                ):

    # k list for one molecule
    k_list = [ round(  i/(num_oscillator-1)+1.0, 4 ) for i in range(num_oscillator) ]
        
    # internal couplings for one molecule
    coupling_mat = [[ 0.0 for i in range(num_oscillator) ] for j in range(num_oscillator) ]

    # fill in internal couplings
    for i in range( num_oscillator ):
        for j in range( i+1, num_oscillator ):
            coupling_mat[i][j] = round( ( k_list[i] - k_list[j] )/ 3.0 , 4 )

    # external couplings for "num_molecule" molecules.
    ext_k_mat = [[ 0.0 for i in range(num_molecule) ] for j in range(num_molecule) ]
    
    for i in range(num_molecule):
        for j in range(num_molecule):
            if i != j:
                ext_k_mat[i][j] = -2.0/ abs(positions[i] - positions[j])**3

                    
    oscillators = [ osys.OscillatorSystem( osys.primitive_oscillator( m = 1.0, k = k_list[i] ) )\
                          for i in range(num_oscillator) ]

    fragments = [ osys.OscillatorSystem( osys.group_of_oscillators(recursive_list = copy.deepcopy( oscillators  ),
                                                                   coupling       = copy.deepcopy( coupling_mat )  ) )
                          for i in range(num_molecule) ]
            
    molecule = osys.OscillatorSystem( osys.group_of_oscillators( recursive_list = fragments,
                                                                 coupling       = ext_k_mat ) )
            
     
    print(molecule.get_k_mat())
    print(molecule.top_level_groups())
    print(molecule.top_level_coupling_mat())
    return molecule

            
def get_lowest_analytical(mol_obj):
    required_states = 5
    k_eff, states_configure, E_analytical, vectors = \
                analytical_solution.analytical_main(
                                                    mol_obj.list_of_masses(),
                                                    mol_obj.get_k_mat(),
                                                    required_states,
                                                    1.0, 1.0, 1.0
                                                    )
    return sorted( E_analytical )[0]
                                                        
                                                        
                                                        
def get_rec_obj( mol_obj, rec_num_states ):
    rec_obj = hamiltonian_class.RecSys_HO_hamiltonian(
                                                      mol_obj.get_k_mat(),
                                                      mol_obj.top_level_groups(),
                                                      mol_obj.top_level_coupling_mat(),
                                                      rec_num_states
                                                      )
    return rec_obj
    

def timeout_handler(signum, frame):
    print("Single Job Max Time Exceeded!")
    raise signal.ItimerError("Single_Job_Time_Out")


def plot_main(
              lanczos_thresh,
              state_per_mol,
              num_oscillator,
              num_molecule,
              inline_positions
              ):
    t_job_start = time()
    thresh = pow(10,lanczos_thresh)
    # Make Total Molecule Obj
    mol_obj = get_mol_obj(
                          num_oscillator,
                          num_molecule,
                          inline_positions      
                          )
    
    # Compute Lowest Analytical Solution
    E_ana = get_lowest_analytical(mol_obj)
    
    # Start iterating Recoupled Solutions by Varying (doubling) "rec_num_states"
    # to achieve above accuracy

    max_time = 7*24*3600 
    
    t_start = time()

    
    output_str = "===================================================\n"

    rec_num_states = [ state_per_mol for i in range(num_molecule) ]
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(max_time)
    try:
        rec_obj =  get_rec_obj( mol_obj, rec_num_states )
        dim     =  rec_obj.get_hamiltonian_dimension()
        size_limit = 128*1024**3//8
        if dim > size_limit:
            output_str += "VECTOR MEMORY REQUEST OVER 128 GB, EXITING..."
            sys.exit(1)
        mat     =  build_matrix(rec_obj)                # imported from test_Lanczos
        vec     =  build_vector(dim)                    # imported from test_Lanczos
        t_rec_start    =  time()
        E_rec          =  experiment(mat,vec,thresh)    # imported from test_Lanczos
        t_rec_end      =  time()
        percent_error  =  abs( E_rec - E_ana ) / E_ana
        
    except signal.ItimerError:
        output_str += "CALCULATION TIME OUT OF " + str(max_time) + " seconds!\n"
        E_rec = -1.0
        t_rec_start    =  time()
        t_rec_end      =  time()
        percent_error  = 100.0
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print(traceback.format_exception(exc_type, exc_value,exc_traceback))
        for item in traceback.format_exception(exc_type, exc_value,exc_traceback):
            output_str += item
    except:
        output_str += "UNKNOWN ERROR TYPE! PLEASE CHECK!\n"
        E_rec = -1.0
        t_rec_start    =  time()
        t_rec_end      =  time()
        percent_error  = 100.0
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print(traceback.format_exception(exc_type, exc_value,exc_traceback))
        for item in traceback.format_exception(exc_type, exc_value,exc_traceback):
            output_str += item
    
    t_exit = time()
    output_str += "===================================================\n"
    output_str += "RESULTS OUTPUT\n"
    output_str += "INPUT INFO\n"
    output_str += "NUM OF STATES  Per MOLECULE = " + str(state_per_mol) + '\n'
    output_str += "NUM OSCILLATOR Per MOLECULE = " + str(num_oscillator) + '\n'
    output_str += "NUM MOLECULE                = " + str(num_molecule) + '\n'

            
    output_str += "POSITIONS Of MOLECULES      = [ "
    for i in inline_positions:
        output_str += str(i) + "  "
    output_str += "]\n"


    output_str += "CALCULATION INFO\n"
    output_str += "Vector   Dimension  = " + str( dim ) + '\n'
    output_str += "REC NUM OF STATES  = " + str(rec_num_states[0]) + '\n'
    output_str += "ANALYTICAL LOWEST EIGENVALUE = " + str(E_ana) + '\n'
    output_str += "RECOUPLE   LOWEST EIGENVALUE = " + str(E_rec) + '\n'
    output_str += "Lanczos Threshold  = %.1e\n" %(thresh)
    output_str += "PERCENT ERROR  = " + str(percent_error) + '\n'
    output_str += "COMPUATATION TIME\n"
    output_str += "Lanczos Extraction Time = "+ str( t_rec_end - t_rec_start ) + " sec\n"
    output_str += "LOOP TOTAL TIME = " + str( t_exit - t_start ) + " sec\n"
    output_str += "JOB  TOTAL TIME = " + str( t_exit - t_job_start ) + " sec\n"
    output_str += "===================================================\n"   
    #print(output_str)
    filename =  "calc_" + str(rec_num_states[0]) + "_states_" + str(num_oscillator) + "_ho_" + str(num_molecule) + "_mol_lan_thresh_"+ str(lanczos_thresh) +".txt"
    print(filename,"Writen to Disk!")
    f = open("/home/yhliu/Qode-clone/Applications/osc-Lanczos/input-output/increase_states_plot/"+ filename,'w')
    f.write(output_str)
    f.close()

    
    



if __name__ == "__main__":
    from multiprocessing import Pool
    from time import ctime
    import os
    
    id = str( os.getpid() )
    job_start_time = ctime()
    
    def wrapper(input_list):
        plot_main( input_list[0], input_list[1], input_list[2], input_list[3], input_list[4])

    lan_thresh = -5
    state_per_mol = [ 2 ]
    num_ho  = [ 2,4,8 ]
    num_mol = [ i for i in range(20,31,2)  ]

    # NUMBER OF THREADS
    num_thrds = 4 
    
    
    input_archieve  = "lanczos threshold     = " + str(lan_thresh) + '\n'
    input_archieve += "states per molecule   = "
    for i in state_per_mol:
        input_archieve += str(i) + "  "
    input_archieve += '\n'
    input_archieve += "number of oscillators = "
    for i in num_ho:
        input_archieve += str(i) + "  "
    input_archieve += '\n'
    input_archieve += "number of molecules   = "
    for i in num_mol:
        input_archieve += str(i) + "  "
    input_archieve += '\n'
    input_archieve += "Job Starting Time = " + job_start_time + '\n'
    input_archieve += "Job pid = " + id + '\n'
    input_archieve += "Num of Threads    = " + str(num_thrds) + '\n'

    archieve_filename = "increase_state_pid_" + id + '_' + job_start_time.replace(" ","_").replace(":","-") + ".inparc"
    f = open("archieve_input/"+ archieve_filename, 'w')
    f.write(input_archieve)
    f.close()



    all_test = []
    
    for n in state_per_mol:
        for i in num_ho:
            for j in num_mol:
                all_test += [[ lan_thresh, n, i, j, [ 10. * k for k in range(j)]  ]]
    

    for line in all_test:
        print(line)
    # pool = Pool(num_thrds)
    pool = Pool(1)
    pool.map( wrapper , all_test )
    pool.terminate()














