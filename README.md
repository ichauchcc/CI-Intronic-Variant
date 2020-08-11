# CI-Intronic-Variant

From vcf file, we get 5 columns of informations. 

Chr - Pos - Ref - Var - strand

This script is to locate any intronic variant which may cause potential splice site change. 

1. get a sequence from reference; use the nucleotide position (+/- 30) (len = 61)
2. If strand == '-'; reverse the sequence
3. Change the Ref to Var in seq No.30 (len-1)/2 
4. Get the sequence 30 +/- 5 to do a sliding window (bin size = 2) check for possible AG/GT
	if AG exist, get the before 20 (+AG) base and after 3 base; count for MaxEnt3
	if GT exist, get the before 3 base and after 6 (+AG) base; count for MaxEnt5
		if MaxEntScore > threshold (5?); count for Ref_MaxEntScore
			if difference > threshold (3?)
				get the sequence, use the nucleotide position (+/- 500) 
				AG backward search 500bp for GT (sliding window size = 2)
				GT forward search 500bp for AG (sliding window size = 2)
				if match, MaxEnt Score it check for threshold
					print the list
