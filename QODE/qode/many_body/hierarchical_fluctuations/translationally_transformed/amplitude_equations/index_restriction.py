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

def get_domain(index, all_operators):
    joint_list = []
    for operator in all_operators:
        if index in operator:
            joint_list += operator
    domain_set = set( joint_list )
    return(domain_set)

def get_dependencies_string(dependencies):
    result  = '\in '
    counter = 0
    if len(dependencies) == 0:
        result += '\infty'
    else:
        for center in dependencies:
            if counter == 0:
                result += 'L_{{{}}}'.format(center)
            else:
                result += '\cap L_{{{}}}'.format(center)
            counter += 1
    return(result)

def get_restrictions(order, all_operators):
    all_operators.sort(key=len)
    past_operators   = []
    all_limits       = []	
    common_indices   = set(order)
    for index in order:
        for operator in all_operators:
            if index in operator:
                domain       = '{}'.format(index)
                index_domain = get_domain(index, all_operators) 
                free_indices = index_domain - common_indices
                dependencies = index_domain & ( set(past_operators) | free_indices )
                #print(index, 'index', dependencies, 'dependencies')
                #result       = get_dependencies_string(dependencies)
                #domain      += result
                #limits_list  = [domain]
                limits_list  = list(dependencies)   # new modification
        all_limits.append(limits_list)
        past_operators.append(index)
        # print( index, result)
    return(all_limits)
	

