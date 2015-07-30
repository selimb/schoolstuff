
from composites.laminate import Laminate
from composites.floatformat import set_options_float
import numpy as np
import scipy.linalg
import sys

def print_shit(a_layup):
    print "Layup : %s " % a_layup
    lam = Laminate(a_layup,1,compute='smart')
    print "V*'s for A: %r" % lam.Vstar_for_A
    print "A* [GPa]:"
    print lam.A/lam.total_ply_thickness
    print "a* [1/TPa]:"
    print lam.a*lam.total_ply_thickness*1000
    print "V*'s for D: %r" % lam.Vstar_for_D
    print "D* : [GPa]:"; Dstar = lam.D/(lam.total_thickness**3/12.0*(1-lam.zc)) 
    print Dstar
    print "d* : [GPa]:"
    print scipy.linalg.inv(Dstar)*1000
    print '\n'

if __name__ == '__main__':
	eval(set_options_float('%f'))
	sys.stdout = open('final_laminates.txt','w')

	print "-----CROSS-PLY------"
	layups = ['0/90s','0/90','0/90/0','0/0/90','0/90/0s','0/0/90s',
			 '90/0/0/90s','90_4/0_4s','90_2/0_2/90_2/0_2s','90/0/90/0/90/0/90/0s'
			 ]
	for layup in layups:
		print_shit(layup)

	print "-----ANGLE-PLY------"
	layups = ['10/-10','10/-10s','60/-60s','-60/60s','45/-45','45/-45s',
				'10/-10/45/-45s','40/-40/60/-60s','40/40/-40s','45/-45/45s']
	for layup in layups:
		print_shit(layup)
	print "-----QUASI-ISO------"
	layups = ['0/-60/60s','0/45/-45/90s','0/36/-36/72/-72s']
	for layup in layups:
		print_shit(layup)


	
	print_shit(layup)
	sys.stdout.close()
