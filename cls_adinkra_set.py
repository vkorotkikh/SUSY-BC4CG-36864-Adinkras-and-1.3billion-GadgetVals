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

	@staticmethod
	def aset_method(tbd_var):
		print(tbd_var, type(tbd_var))

		# pass

	def __init__(self, dim, nodes, adink_list):

		self.dim 		= dim
		self.nodes		= nodes

		self.adinkra_list = adink_list


	def ret_adinkras(self):		# class instance method
		# what will this do?
		return self.adinkra_list

	def get_len_adinkra(self):
		return len(self.adinkra_list)
