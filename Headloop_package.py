#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package for designing headloop primers to use in headloop suppression PCR to 
suppress amplification of a known haplotype.
Copyright (C) July 2020, Gareth T. Powell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


Function designs tags from guide sequence provided by the user, with frameshifting 
to minimise Tm differences between the tags and the PCR primers, then adds the 'best' 
tags to the correct PCR primer provided by the user, depending upon strand orientation 
of the guide sequence relative to the PCR primers e.g. if the guide sequence in on the
same strand as the sense primer, a reverse complement tag is added to the sense 
primer and an offset tag is added to the antisense primer.

For further details of headloop suppression PCR, see the preprint at bioRxiv:

    Kroll, F., Powell, G.T. et al (2021) eLife 10:e59683, doi: 10.7554/eLife.59683
    https://doi.org/10.7554/eLife.59683

This package uses 'melting' from Erik Clarke [https://github.com/eclarke/melt] to 
calculate Tm of primers and headloop tags, and Seq from Biopython 
[https://biopython.org/wiki/Seq] for reverse complementation.

To use this package, the user needs to provide four variables:
    sense_oligo - string containing the sense primer
    antisense_oligo - string containing the antisense primer
    guide_context - string containing guide sequence and >= 15 bp forward context
    orientation - is the guide in the same strand as the 'sense' primer or 'antisense' 
                  primer?

Example (tbx16_AA):
  design('CTGGTCCAGTGCGTTATTGG', 'AGCCAAATGCTTCTTGCTCTTTT', 
           'CTACAGGACGTACCTGCACCCGGATTCACCAGCGCCCG', 'antisense')
    
Returns:
    Two primers as SeqRecord objects, with comments on Tm matching in the description
    
    CCTGCACCCGGATTCACCAGCTGGTCCAGTGCGTTATTGG
    WARNING: Could not optimise sense headloop tag (Tm difference > 3°C) 
    
    GGTGCAGGTACGTCCTGTAGAGCCAAATGCTTCTTGCTCTTTT
    Tm difference < 3°C
    
    (Gives a warning flag if the Tm difference between the headloop tag and the 
     base primer is calculated to be > 3C)
    
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import melting #library from Erik Clarke's GitHub repository [https://github.com/eclarke/melt]

def design(sense_oligo, antisense_oligo, guide_context, orientation):

    #identify strand for Tm comparison and tagging   
    orient_error = 'Orientation not correctly specified.\n'
    
    if orientation == 'sense':
        i = 0
        j = 1
    elif orientation == 'antisense':
        i = 1
        j = 0
    else:
        print(orient_error)
        return
    
    #convert sequences to Seq objects for manipulation    
    sense = Seq(sense_oligo)
    antisense = Seq(antisense_oligo)
    context = Seq(guide_context)

    primers = (sense, antisense)
    
    #check guide sequence + context is big enough for design
    size_error = 'Guide and context is not big enough for design'
    
    if len(context) < 35:
        print(size_error)
        return
    
    #calculate Tm for primers
    MT = []

    for primer in primers:
        temptm = melting.temp(primer, DNA_c = 1000, Mg_c = 1.5, dNTPs_c = 0.2)
        MT.append(temptm)
    
    #design headloop tags
    frame = 0
    poss_rc = []    #initialise lists
    poss_offset = []
    master_rc = []  #initialise master lists
    master_offset = []
    
    while frame < 3:
        
        #create list of possible reverse complement tags
        temp_rc = context[frame:(frame + 20)]    #extract sequence
        rc = temp_rc.reverse_complement()    #reverse complement temp_rc
        rc_temp_tm = melting.temp(rc, DNA_c = 1000, Mg_c = 1.5, dNTPs_c = 0.2)     #calculate Tm melting method
        rc_temp_diff = abs(MT[i] - rc_temp_tm)  #calculate Tm difference with sense primer
        if rc_temp_diff < 3:
            rc_flag = 0
        else:
            rc_flag = 1
        poss_rc = (rc, rc_temp_tm, rc_temp_diff, frame, rc_flag)    #create list
        master_rc.append(poss_rc) #append to master list

        #create list of possible offset tags
        offset = context[(frame + 12):(frame + 32)]     #extract offset sequence
        off_temp_tm = melting.temp(offset, DNA_c = 1000, Mg_c = 1.5, dNTPs_c = 0.2)     #calculate Tm melting method
        off_temp_diff = abs(MT[j] - off_temp_tm)  #calculate Tm difference with antisense primer
        if off_temp_diff <3:
            off_flag = 0
        else:
            off_flag = 1
        poss_offset = (offset, off_temp_tm, off_temp_diff, frame, off_flag)  #create list
        master_offset.append(poss_offset) #append to master list
        frame = frame + 1 #frameshift by 1
    
    #sort lists of optional tags by Tm difference with primer and, in tiebreaks, position
    master_rc.sort(key = lambda x: (x[2], x[3]))
    master_offset.sort(key = lambda x: (x[2], -x[3]))
    
    #concatenate headloop tags with correct primer, depending on strand orientation
    if orientation == 'sense':
        temp_sense_hl = master_rc[0][0] + sense #create headloop sense primer sequence
        temp_antisense_hl = master_offset[0][0] + antisense #create headloop antisense primer sequence
        
        #formats seq records
        sense_hl = SeqRecord(temp_sense_hl, id = 'Sense HL')
        if master_rc[0][4] == 1:
            sense_hl.description = 'WARNING: Could not optimise sense headloop tag (Tm difference > 3\u00b0C)'
        else:
            sense_hl.description = 'Tm difference < 3\u00b0C'

        antisense_hl = SeqRecord(temp_antisense_hl, id = 'Antisense HL')
        if master_offset[0][4] == 1:
            antisense_hl.description = 'WARNING: Could not optimise antisense headloop tag (Tm difference > 3\u00b0C)'
        else:
            antisense_hl.description = 'Tm difference < 3\u00b0C'

            
    else:
        temp_antisense_hl = master_rc[0][0] + antisense #create headloop antisense primer sequence
        temp_sense_hl = master_offset[0][0] + sense

        #formats seq records
        sense_hl = SeqRecord(temp_sense_hl, id = 'Sense HL')
        if master_offset[0][4] == 1:
            sense_hl.description = 'WARNING: Could not optimise sense headloop tag (Tm difference > 3\u00b0C)'
        else:
            sense_hl.description = 'Tm difference < 3\u00b0C'
            
        antisense_hl = SeqRecord(temp_antisense_hl, id = 'Antisense HL')
        if master_rc[0][4] == 1:
            antisense_hl.description = 'WARNING: Could not optimise antisense headloop tag (Tm difference > 3\u00b0C)'
        else:
            antisense_hl.description = 'Tm difference < 3\u00b0C'
            
    return (sense_hl, antisense_hl)
