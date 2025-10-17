# dist_calc
R script to compute intra- and inter-MOTU distances after SWARM
Only MOTUs with more than one ESV are used
For inter-MOTU distances it uses the representative sequence of each MOTU

The input file must be a .csv with all ESVs in rows. Columns can contain any information (typically values of abundance per samples), but MUST include a column "MOTU" with the code of the MOTU to which
each ESV is assigned, and a column "seq" with the sequence of the ESV (all sequences should be aligned). 
ESVs must be ordered by abundance, so the first ESV of each MOTU has the representative sequence
