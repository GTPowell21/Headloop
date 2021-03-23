#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script to demonstrate use of headloop package to design primers for headloop suppression PCR.
Copyright (C) July 2020, Gareth T. Powell
"""

from headloop.designer import design

##Kroll et al 2021 eLife 10:e59683, https://doi.org/10.7554/eLife.59683
##Example from the paper: tbx16_AA

#Forward primer
tbx16_AA_F = 'AGGTTATTTGCTGTCATGGCTTTG'

#Reverse primer
tbx16_AA_R = 'ACTTTCACATCATTCCACTGG'

#Guide sequence and PAM, plus 15 bp downstream sequence
tbx16_AA_guide = 'ACCATCATGTGCTGGACGTCCGGATTGATGGAGCG'

#Orientation of guide with respect to forward primer: same strand = sense, opposite strand = 'antisense'
tbx16_AA_guide_orientation = 'antisense'

#Run the function
output = design(tbx16_AA_F, tbx16_AA_R, tbx16_AA_guide, tbx16_AA_guide_orientation)

#Output of the function: note that 'design' returns two primers as SeqRecord objects, with comments on how closely matched
#the predicated annealing temperatures are. Only one of those primers should be used, to be determined by testing.
#In this individual case, the antisense headloop primer was used for HL-PCR.
print(' Sense headloop primer:', output[0].seq, '\n', output[0].description, '\n', 'Antisense headloop primer:', output[1].seq, '\n', output[1].description)