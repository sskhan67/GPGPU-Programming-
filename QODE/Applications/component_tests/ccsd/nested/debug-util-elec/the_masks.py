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
from qode.many_body.nested_operator import mask

K   = [True, False, False, False]
E   = [False, K, False, False ]
F   = [False, False, K, False ]
D   = [False, False, False, K ]
EE  = [False, E, False, False ]
EF  = [False, F, False, False ]
ED  = [False, D, False, False ]
FF  = [False, False, F, False ]
FD  = [False, False, D, False ]
DD  = [False, False, False, D ]
EED = [False, ED, False, False ]
EFF = [False, FF, False, False ]
EFD = [False, FD, False, False ]

Emask   = mask(E)
Fmask   = mask(F)
Dmask   = mask(D)
EEmask  = mask(EE)
EFmask  = mask(EF)
EDmask  = mask(ED)
FFmask  = mask(FF)
FDmask  = mask(FD)
DDmask  = mask(DD)
EEDmask = mask(EED)
EFFmask = mask(EFF)
EFDmask = mask(EFD)

Xmasks = {\
"F"   : Fmask,
"D"   : Dmask,
"EF"  : EFmask,
"ED"  : EDmask,
"FF"  : FFmask,
"FD"  : FDmask,
"DD"  : DDmask,
"EED" : EEDmask,
"EFF" : EFFmask,
"EFD" : EFDmask
}

Tmasks = {\
"E"   : Emask,
"EE"  : EEmask
}

all_masks = { "K":mask(K) }
all_masks.update(Xmasks)
all_masks.update(Tmasks)
