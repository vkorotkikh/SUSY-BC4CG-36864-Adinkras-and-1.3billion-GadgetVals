# ******************************************************************************
# Name:    Caculate all BC4-based Adinkras
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    September 2017
# Version: N/A
#
# Description: The code calculates all all unique 36,864 ordered BC4-based
# adinkras with four colors, four open-nodes and four closed nodes and then
# calculates the Vij holoraumy matrices and the Gadget values
#
# ******************************************************************************

def main():
	pass



class AdinkraSet():

	@classmethod
	def aset_classmethod(cls):
		print("Class Name: %s" % cls.__name__)

	def __init__(self, dim, nodes, adink_list):

		self.dim 		= dim
		self.nodes		= nodes

		self.adinkra_list = adink_list


	def rand_inst_method(self):		# class instance method
		# what will this do?
		return self.dim

	def len_adinkralist(self):
		return len(self.adinkra_list)