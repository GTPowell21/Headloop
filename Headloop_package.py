#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package for designing headloop primers to use in headloop suppression PCR to 
suppress amplification of a known haplotype.
Copyright (C) July 2020, Gareth T. Powell

Function designs tags from guide sequence provided by the user, with frameshifting 
to minimise Tm differences between the tags and the PCR primers, then adds the 'best' 
tags to the correct PCR primer provided by the user, depending upon strand orientation 
of the guide sequence relative to the PCR primers e.g. if the guide sequence in on the
same strand as the sense primer, a reverse complement tag is added to the sense 
primer and an offset tag is added to the antisense primer.

For further details of headloop suppression PCR, see the preprint at bioRxiv:

                  https://doi.org/10.1101/2020.06.04.133462

This package uses 'melting' from Erik Clarke [https://github.com/eclarke/melt] to 
calculate Tm of primers and headloop tags, and Seq from Biopython 
[https://biopython.org/wiki/Seq] for reverse complementation.

To use this package, the user needs to provide four variables:
    sense_oligo - string containing the sense primer
    antisense_oligo - string containing the antisense primer
    guide_context - string containing guide sequence and >= 15 bp forward context
    orientation - is the guide in the same strand as the 'sense' primer or 'antisense' 
                  primer?

Example:
    headloop.design('ATCG', 'CGAT', 'AGCTCTGT...TGG', 'sense')
    
Returns:
    'HL+ATCG', 'HL+CGAT', sense_warning?, antisense_warning?
    
    (Gives a warning flag if the Tm difference between the headloop tag and the 
     base primer is calculated to be > 3C)
    
"""

from Bio import Seq
import melting #library from Erik Clarke's GitHub repository [https://github.com/eclarke/melt]

def design(sense_oligo, antisense_oligo, guide_context, orientation):

    #identify strand for Tm comparison and tagging   
    orient_error = 'Orientation not correctly specified'
    
    if orientation == 'sense':
        switch = -1
    elif orientation == 'antisense':
        switch = 0
    else:
        return orient_error
    
    #convert sequences to Seq objects for manipulation    
    sense = Seq(sense_oligo)
    antisense = Seq(antisense_oligo)
    context = Seq(guide_context)

    primers = (sense, antisense)
    
    #check guide sequence + context is big enough for design
    size_error = 'Guide and context is not big enough for design'
    
    if (len(context.seq) + 1) < 35:
        return size_error
    
    #calculate Tm for primers
    MT = []

    for primer in primers:
        temptm = melting.temp(primer.seq, DNA_c = 1000, Mg_c = 1.5, dNTPs_c = 0.2)
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
        rc_temp_tm = melting.temp(rc.seq, DNA_c = 1000, Mg_c = 1.5, dNTPs_c = 0.2)     #calculate Tm melting method
        rc_temp_diff = abs(MT[(switch + 1)] - rc_temp_tm)  #calculate Tm difference with sense primer
        if rc_temp_diff < 3:
            rc_flag = 0
        else:
            rc_flag = 1
        poss_rc = (rc, rc_temp_tm, rc_temp_diff, frame, rc_flag)    #create list
        master_rc.append(poss_rc) #append to master list

        #create list of possible offset tags
        offset = context[(frame + 12):(frame + 32)]     #extract offset sequence
        off_temp_tm = melting.temp(offset.seq, DNA_c = 1000, Mg_c = 1.5, dNTPs_c = 0.2)     #calculate Tm melting method
        off_temp_diff = abs(MT[(switch + 2)] - off_temp_tm)  #calculate Tm difference with antisense primer
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
        sense_hl = master_rc[0][0] + sense #create headloop sense primer sequence
        antisense_hl = master_offset[0][0] + antisense #create headloop antisense primer sequence
        
        #formats seq records
        sense_hl.id = 'Sense HL'
        if master_rc[0][4] == 1:
            sense_hl.description = 'WARNING: Could not optimise sense headloop tag (Tm difference > 3\u00b0C)'

        antisense_hl.id = 'Antisense HL'
        if master_offset[0][4] == 1:
            antisense_hl.description = 'WARNING: Could not optimise antisense headloop tag (Tm difference > 3\u00b0C)'
            
    else:
        antisense_hl = master_rc[0][0] + antisense #create headloop antisense primer sequence
        sense_hl = master_offset[0][0] + sense

        #formats seq records
        sense_hl.id = 'Sense HL'
        if master_offset[0][4] == 1:
            sense_hl.description = 'WARNING: Could not optimise sense headloop tag (Tm difference > 3\u00b0C)'

        antisense_hl.id = 'Antisense HL'
        if master_rc[0][4] == 1:
            antisense_hl.description = 'WARNING: Could not optimise antisense headloop tag (Tm difference > 3\u00b0C)'

    return (sense_hl, antisense_hl)

