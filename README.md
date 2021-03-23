# Headloop
Package for designing headloop primers to use in headloop suppression PCR to 
suppress amplification of a known haplotype.
Copyright (C) July 2020, Gareth T. Powell

## License
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

## Description
Function designs tags from guide sequence provided by the user, with frameshifting 
to minimise Tm differences between the tags and the PCR primers, then adds the 'best' 
tags to the correct PCR primer provided by the user, depending upon strand orientation 
of the guide sequence relative to the PCR primers e.g. if the guide sequence in on the
same strand as the sense primer, a reverse complement tag is added to the sense 
primer and an offset tag is added to the antisense primer.

For further details of headloop suppression PCR, see the paper at eLife:
    
    Kroll, F., Powell, G.T. et al (2021) eLife 10:e59683, doi: 10.7554/eLife.59683
    https://doi.org/10.7554/eLife.59683

## Installation

    pip install headloop

This package uses 'melting' from Erik Clarke [https://github.com/eclarke/melt] to 
calculate Tm of primers and headloop tags, and Seq & SeqRecord from Biopython 
[https://biopython.org/wiki/Seq] for reverse complementation and object input/output.

## Usage
To use this package, the user needs to provide four variables:
    
    sense_oligo     #string containing the sense primer
    antisense_oligo #string containing the antisense primer
    guide_context   #string containing guide sequence and >= 15 bp forward context
    orientation     #is the guide in the same strand as the 'sense' primer or 'antisense' primer?

Example (tbx16_AA):

    from headloop.designer import design
    
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
