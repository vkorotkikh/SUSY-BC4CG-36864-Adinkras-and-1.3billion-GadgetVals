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
	def adinkra_clsmethod(cls):
		print("Class Name: %s" % cls.__name__)

	def __init__(self, dim, nodes):

		self.dim = dim
		self.nodes = nodes


	def rand_inst_method(self):		# class instance method
		# what will this do?
		return self.dim
