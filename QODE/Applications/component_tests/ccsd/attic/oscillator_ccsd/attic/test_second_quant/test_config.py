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
# from Applications.general_ci import config, state, generate_config
from Applications.general_ci import config

# One Vacuum Config
config_vaccum = [ False, False, False, False]


#
# Testing Creations
#
print("Creations")
# Initialize a configuration object
config_1    = config.configuration(config_vaccum, coeff=1.0)
config_1.print_info()

# Create at orbital 0
config_2 = config.create(0, config_1)
config_2.print_info()

# Create at orbital 1
config_3 = config.create(1, config_2)
config_3.print_info()

# Create at orbital 3
config_4 = config.create(3, config_3)
config_4.print_info()

# Create at orbital 1 again to get None.
config_5 = config.create(1, config_4)
config_5.print_info()

# Create at a None config one more time, we should get a None.
config_6 = config.create(0, config_5)
config_6.print_info()



#
# Starting a new thread for testing annihilations
#
print("Annihilations")

# Initialize a configuration object
config_7    = config.configuration(config_vaccum, coeff=1.0)
config_7.print_info()

# Create at orbital 0
config_8 = config.create(0, config_7)
config_8.print_info()

# Create at orbital 1
config_9 = config.create(1, config_8)
config_9.print_info()

# Create at orbital 3
config_10 = config.create(3, config_9)
config_10.print_info()

# Annihilate orbital 0
config_11 = config.annihilate(1, config_10)
config_11.print_info()

# Annihilate orbital 3
config_12 = config.annihilate(0, config_11)
config_12.print_info()

# Annihilate orbital 0 again, we should get a None
config_13 = config.annihilate(0, config_12)
config_13.print_info()

# Annihilate at None to get another None
config_14 = config.annihilate(0, config_13)
config_14.print_info()







