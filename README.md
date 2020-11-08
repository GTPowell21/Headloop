# Headloop
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

Kroll et al. 2020 [https://doi.org/10.1101/2020.06.04.133462]

This package uses 'melting' from Erik Clarke [https://github.com/eclarke/melt] to 
calculate Tm of primers and headloop tags, and Seq from Biopython 
[https://biopython.org/wiki/Seq] for reverse complementation.

To use this package, the user needs to provide four variables:
    sense_oligo - string containing the sense primer;
    antisense_oligo - string containing the antisense primer;
    guide_context - string containing guide sequence and >= 15 bp forward context;
    orientation - is the guide in the same strand as the 'sense' primer or 'antisense' 
                  primer?

Example:
    headloop.design('ATCG', 'CGAT', 'AGCTCTGT...TGG', 'sense')
    
Returns:
    'HL+ATCG', 'HL+CGAT', sense_warning?, antisense_warning?
    
   (Gives a warning flag if the Tm difference between the headloop tag and the 
    base primer is calculated to be > 3C)
