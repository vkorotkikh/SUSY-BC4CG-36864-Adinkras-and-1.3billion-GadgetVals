# ******************************************************************************
# Name:    Caculate All Possible BC4-based Adinkras
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    December 2016
# Version: 1.3A
#
# Description: The code calculates all all unique 36,864 ordered BC4-based
# adinkras with four colors, four open-nodes and four closed nodes and then
# calculates the Vij holoraumy matrices and the Gadget values
#
# ******************************************************************************

# ******************************************************************************
# Library Imports
import sys, time, math, itertools
import numpy as np
import numpy.matlib


# ******************************************************************************
# Function Imports
import adinkra_nxn_constructor
import cls_adinkra_set

# ******************************************************************************
# Main() function.
def main():
	# Just making standard 4 4 for now. Quick OO test build
	adinkra_list	= []
	adinkra_list	= adinkra_nxn_constructor.makeall_adinkras(4,4)
	print("Length adinkra_list: ", len(adinkra_list))

	cls_adinkra_set.AdinkraSet.aset_classmethod()
	if len(adinkra_list) > 1:
		# NewAdink = cls_adinkra_set.AdinkraSet(4,4, adinkra_list)
		NewAdink = cls_adinkra_set.AdinkraSet(4,4,adinkra_list)
		print("Len Adinkra Class list: ", NewAdink.get_len_adinkra())


if __name__ == "__main__":
	start_time = time.time()
	main()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
