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
            
     
    #print(molecule.get_k_mat())
    #print(molecule.top_level_groups())
    #print(molecule.top_level_coupling_mat())
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

def prep_for_search( distance ):
    # It is a BINARY SEARCH, divided by 2.
    div_result = distance // 2
    num_search = 0
    shifts     = []
    while div_result > 0:
        #print( "Result =",div_result )
        num_search += 1
        shifts += [ div_result ]
        div_result = div_result // 2
    print("NUM OF MAX SEARCH =",num_search)
    print("SHIFTS OF  SEARCH =",shifts)
    return num_search, shifts

            
def plot_main(
              accuracy_required,
              lanczos_thresh,
              num_oscillator,
              num_molecule,
              inline_positions
              ):
    # Calculate Accuracy Range:
    t_job_start = time()
    accuracy_range = [ accuracy_required -0.5, accuracy_required + 0.5 ]
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
    rec_num_states = [ 2 for i in range(num_molecule) ]
    run = True

    max_time = 3600*200
    
    rec_state_history = []
    output_str  = ""
    num_iter = 0
    high_accu_limit = pow(10, accuracy_range[0])
    low_accu_limit  = pow(10, accuracy_range[1])
    t_start = time()
    while run:
        num_iter += 1
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(max_time)
        try:
            rec_obj =  get_rec_obj( mol_obj, rec_num_states )
            dim     =  rec_obj.get_hamiltonian_dimension()
            mat     =  build_matrix(rec_obj)                # imported from simple_fast
            vec     =  build_vector(dim)                    # imported from simple_fast
            t_rec_start    =  time()
            E_rec          =  experiment(mat,vec,thresh)      # imported from simple_fast
            #output_str = "CALCULATION TIMEOUT FOR 3 HOURS\n\n" + output_str
            t_rec_end      =  time()
            percent_error  =  abs( E_rec - E_ana ) / E_ana
            rec_state_history += [ rec_num_states[0] ]
            
            print("===================================================")
            print("Checking Recouple Solution Accuracy...")
            if high_accu_limit <  percent_error <  low_accu_limit:
                print("Iteration",num_iter,"Accuracy In Range, Normal Termination")
                output_str += "Iteration "+ str(num_iter) + " Accuracy In Range, Normal Termination\n"
                run = False
            elif percent_error > low_accu_limit:
                if rec_state_history[-1] == max(rec_state_history):
                    print("Iteration",num_iter,"Accuracy Too Low, Doubling States")
                    output_str += "Iteration "+ str(num_iter) + " Accuracy Too Low, Doubling States\n"
                    for i in range(num_molecule):
                        rec_num_states[i] *= 2
                else:
                    # Continue Binary Search
                    if current_search < num_search:
                        for i in range(num_molecule):
                            rec_num_states[i] += shifts[current_search]
                        current_search += 1
                        print("Iteration",num_iter,"Accuracy Too Low, Continue Binary Search")
                        output_str += "Iteration "+ str(num_iter) + " Accuracy Too Low, Continue Binary Search\n"
                        
                    else:
                        # Highest Number has a too Low accuracy, Go back to any great number than current num of states.
                        for i in range(num_molecule):
                            rec_num_states[i] = search_range[-1]
                        rec_state_history += [ rec_num_states[0] ]
                        rec_obj =  get_rec_obj( mol_obj, rec_num_states )
                        dim     =  rec_obj.get_hamiltonian_dimension()
                        mat     =  build_matrix(rec_obj)                # imported from simple_fast
                        vec     =  build_vector(dim)                    # imported from simple_fast
                        t_rec_start    =  time()
                        E_rec          =  experiment(mat,vec,thresh)      # imported from simple_fast
                        #output_str = "CALCULATION TIMEOUT FOR 3 HOURS\n\n" + output_str
                        t_rec_end      =  time()
                        percent_error  =  abs( E_rec - E_ana ) / E_ana
                        # Terminate the Loop 
                        run = False
                        print("Iteration",num_iter,"Accuracy Too Low, Go Back to Search Range High Limit")
                        output_str += "Iteration "+ str(num_iter) + " Accuracy Too Low, Go Back to Search Range High Limit"
                            
            elif percent_error < high_accu_limit and rec_num_states[0] > 2:
                if rec_state_history[-1] == max(rec_state_history):
                    print("Iteration",num_iter,"Accuracy Too High, Reducing States")
                    output_str += "Iteration "+ str(num_iter) + " Accuracy Too High, Reducing States\n"
                    # Start A Binary Search
                    distance = rec_state_history[-1] - rec_state_history[-2]
                    print("Search Range =",rec_state_history[-2:], "DISTACNE =", distance)
                    search_range = [ rec_state_history[-2], rec_state_history[-1] ]
                    num_search, shifts = prep_for_search( distance )
                    current_search = 0
                    for i in range(num_molecule):
                        rec_num_states[i] -= shifts[current_search]
                    current_search += 1
                else:
                    # Continue Binary Search.
                    print("Iteration",num_iter,"Accuracy Too High, Reducing States")
                    output_str += "Iteration "+ str(num_iter) + " Accuracy Too High, Reducing States\n"
                    if current_search < num_search:
                        for i in range(num_molecule):
                            rec_num_states[i] -= shifts[current_search]
                        current_search += 1
                    else:
                        # Lowest Number hsa a too high accuracy, just keep it as final.
                        print("Iteration ",num_iter," Accuracy Too High for the Least possible state number... Exiting...")
                        output_str += "Iteration "+str(num_iter)+" Accuracy Too High for the Least Number in Binary Search... Num of States Found\n"
                        run = False


            else:
                print("Iteration",num_iter,"Accuracy Too High for Even Two States, Exiting...")
                output_str += "Iteration "+ str(num_iter) + " Accuracy Too High for Even Two States, Exiting...\n"
                run = False
        except signal.ItimerError:
            print("CALCULATION TIME OUT OF %d SECONDS" %(max_time))
            output_str += "CALCULATION TIME OUT of %d SECONDS\n" %(max_time)
            print("Unexpected error:", sys.exc_info())
            t_rec_start    =  time()
            t_rec_end      =  time()
            E_rec = -1.0
            percent_error = -100.0
            run = False
        
        except:
            print("JOB ERROR TERMINATION")
            output_str += "JOB ERROR TERMINATION\n"
            t_rec_start    =  time()
            t_rec_end      =  time()
            E_rec = -1.0
            percent_error = -100.0
            run = False
            exc_type, exc_value, exc_traceback = sys.exc_info()
            #print("Unexpected error:", sys.exc_info())
            #print(traceback.print_exc())
            print(traceback.format_exception(exc_type, exc_value,exc_traceback))
            for item in traceback.format_exception(exc_type, exc_value,exc_traceback):
                output_str += item
        print("===================================================")

    
    t_exit = time()
    output_str += "===================================================\n"
    output_str += "OUTPUT RESULTS\n"
    output_str += "NUM ITERATIONS = " + str(num_iter) + '\n'
    output_str += "NUM MOLECULE   = " + str(num_molecule) + '\n'
    output_str += "NUM OSCILLATOR = " + str(num_oscillator) + '\n'
    output_str += "POSITIONS      = [ "
    for i in inline_positions:
        output_str += str(i) + "  "
    output_str += "]\n"
    output_str += "States History = [ "
    for  i in rec_state_history:
        output_str += str(i) + "  "
    output_str += "]\n"
    output_str += "ACCURACY REQUAIRED = 1.0e" + str(accuracy_required) + '\n'
    output_str += "Lanczos Threshold  = %.1e\n" %(thresh)
    output_str += "Vector   Dimesion  = " + str( dim ) + '\n'
    output_str += "Actual Rec States  = " + str(rec_num_states[0]) + '\n'
    output_str += "ANALYTICAL LOWEST EIGENVALUE = " + str(E_ana) + '\n'
    output_str += "RECOUPLE   LOWEST EIGENVALUE = " + str(E_rec) + '\n'
    output_str += "ACCURACY RANGE = %.6e ~ %.6e\n" %(  low_accu_limit  , high_accu_limit )
    output_str += "PERCENT ERROR  = " + str(percent_error) + '\n'
    output_str += "Lanczos Extraction Time = "+ str( t_rec_end - t_rec_start ) + " sec\n"
    output_str += "LOOP TOTAL TIME = " + str( t_exit - t_start ) + " sec\n"
    output_str += "JOB  TOTAL TIME = " + str( t_exit - t_job_start ) + " sec\n"
    output_str += "===================================================\n"   
    #print(output_str)
    filename = "accu_" + str(accuracy_required) + "_calc_" + str(num_molecule) + "_mol_" + str(num_oscillator) + "_ho.txt"
    print(filename,"Writen to Disk!")
    f = open("/home/yhliu/Qode-clone/Applications/osc-Lanczos/input-output/fix_accuray_plots/"+ filename,'w')
    f.write(output_str)
    f.close()

    
    



if __name__ == "__main__":
    from multiprocessing import Pool
    
    def wrapper(input_list):
        plot_main( input_list[0],input_list[1],input_list[2],input_list[3],input_list[4] )

    accu_list =  [-9]
    thresh_list = [-5]
    num_ho  = [ 2 ]
    num_mol = [ i for i in range(7, 25) ]

    
    all_test = []
    for n in range(len(accu_list)):
        for i in num_ho:
            for j in num_mol:
                all_test += [[ accu_list[n] , thresh_list[n] , i, j, [ 10. * k for k in range(j)]  ]]
    for line in all_test:
        print(line)
    pool = Pool(1)
    
    pool.map( wrapper , all_test )















