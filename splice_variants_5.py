def genome_ref(chr_pos,ref,strand):
    if strand == '+':
        return Seq(str(ref[chr_pos][0:]))
    if strand == '-':
        return Seq(str(ref[chr_pos][0:])).reverse_complement()
      
# =============================================================================
#     alt_map = {'ins':'0'}
#     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
# def reverse_comp(seq):    
#     for k,v in alt_map.iteritems():
#         seq = seq.replace(k,v)
#     bases = list(seq) 
#     bases = reversed([complement.get(base,base) for base in bases])
#     bases = ''.join(bases)
#     for k,v in alt_map.iteritems():
#         bases = bases.replace(v,k)
#     return bases
# =============================================================================


############################################################################
#Check splice donor or acceptor
def splice_seq(var_pos,var_b,ref):
    if var_b.upper() == 'A':
        if ref[var_pos+1].upper() == 'G':
            return 'Acceptor', str(ref[var_pos-18:var_pos] + var_b + ref[var_pos+1:var_pos+5]), str(ref[var_pos-18:var_pos+5])
    if  var_b.upper() == 'T':
        if ref[var_pos-1].upper() == 'G':
            return 'Donor', str(ref[var_pos-4:var_pos] + var_b + ref[var_pos+1:var_pos+5]), str(ref[var_pos-4:var_pos+5])
    if var_b.upper() == 'G':
        if ref[var_pos-1].upper() == 'A':
            return 'Acceptor', str(ref[var_pos-19:var_pos] + var_b + ref[var_pos+1:var_pos+4]), str(ref[var_pos-19:var_pos+4])
        if ref[var_pos+1].upper() == 'T':
            return 'Donor', str(ref[var_pos-3:var_pos] + var_b + ref[var_pos+1:var_pos+6]), str(ref[var_pos-3:var_pos+6])


#Input variant to different function according to reference genome direction
def splice_var(var_pos,var_b,strand,ref):
    if strand == '+':
        return splice_seq(var_pos-1,var_b, ref)
    if strand == '-':
        return splice_seq(len(ref)-var_pos,var_b,ref)


############################################################################      
#If splice variant is an ACCEPTOR, find adjacent DONOR
def donor_adj(width, var_pos, var_b,ref,strand):
    donor_adj_dict = {}
    if var_b.upper() == 'G':
        search_seq = ref[var_pos-width-1:var_pos-1]
        for base in range(19,len(search_seq)-4):
            if search_seq[base].upper() == 'A':
                if search_seq[base+1].upper() == 'G':
                    if strand == '+':
                        donor_adj_dict[search_seq[base-18:base+5]] = {'Start':var_pos-len(search_seq)+base-19, 'End':var_pos-len(search_seq)+base+3}
                    if strand == '-':
                        donor_adj_dict[search_seq[base-18:base+5]] = {'Start':len(ref)-var_pos+len(search_seq)-base+18,'End':len(ref)-var_pos+len(search_seq)-base-4}
    if var_b.upper() == 'T':
        search_seq = ref[var_pos-width-2:var_pos-2]
        for base in range(19,len(search_seq)-4):
            if search_seq[base].upper() == 'A':
                if search_seq[base+1].upper() == 'G':
                    if strand == '+':
                        donor_adj_dict[search_seq[base-18:base+5]] = {'Start':var_pos-len(search_seq)+base-19, 'End':var_pos-len(search_seq)+base+3}
                    if strand == '-':
                        donor_adj_dict[search_seq[base-18:base+5]] = {'Start':len(ref)-var_pos+len(search_seq)-base+18,'End':len(ref)-var_pos+len(search_seq)-base-4}
    return donor_adj_dict
                    

#If splice variant is a DONOR, find adjacent ACCEPTOR
def acceptor_adj(width, var_pos, var_b,ref,strand):
    acceptor_adj_dict = {}
    if var_b.upper() == 'A':
        search_seq = ref[var_pos+2:var_pos+width+2]
        for base in range(4,len(search_seq)-5):
            if search_seq[base].upper() == 'G':
                if search_seq[base+1].upper() == 'T':
                    if strand == '+':
                        acceptor_adj_dict[search_seq[base-3:base+6]] = {'Start':var_pos-len(search_seq)+base-4, 'End':var_pos-len(search_seq)+base+4}
                    if strand == '-':
                        acceptor_adj_dict[search_seq[base-3:base+6]] = {'Start':len(ref)-var_pos+len(search_seq)-base+2,'End':len(ref)-var_pos+len(search_seq)-base-6}
    if var_b.upper() == 'G':
        search_seq = ref[var_pos+1:var_pos+width+1]
        for base in range(4,len(search_seq)-5):
            if search_seq[base].upper() == 'G':
                if search_seq[base+1].upper() == 'T':
                    if strand == '+':
                        acceptor_adj_dict[search_seq[base-3:base+6]] = {'Start':var_pos-len(search_seq)+base-4, 'End':var_pos-len(search_seq)+base+4}
                    if strand == '-':
                        acceptor_adj_dict[search_seq[base-3:base+6]] = {'Start':len(ref)-var_pos+len(search_seq)-base+2,'End':len(ref)-var_pos+len(search_seq)-base-6}
    return acceptor_adj_dict


#Input variant to different function according to splice nature
def adj_var_seq(var_pos,var_b,strand,ref,nature,width):
    if strand == '+':
        if nature.upper() == 'DONOR':
            return donor_adj(width, var_pos, var_b,ref,strand)
        if nature.upper() == 'ACCEPTOR':
            return acceptor_adj(width, var_pos-1, var_b,ref,strand)
    if strand == '-':
        if nature.upper() == 'DONOR':
            return donor_adj(width, len(ref)-var_pos, var_b,ref,strand)
        if nature.upper() == 'ACCEPTOR':
            return acceptor_adj(width, len(ref)-var_pos, var_b,ref,strand)
      

############################################################################
def splice_search(win_seq,var_pos,strand,ref,win_width):
    splice_dict = {}
    if strand == '+':
        for base in range(25,len(win_seq)-26):
            if win_seq[base].upper() == 'G':
                if win_seq[base-1].upper() == 'A':
                    potential_seq = str(win_seq[base-19:base+4])
                    maxentscore_pot = maxent_fast.score3(potential_seq, matrix=matrix3)
                    if maxentscore_pot > 2:
                        splice_dict[potential_seq] = {'Nature':'Acceptor','Start':var_pos-30+base-19,'End':var_pos-30+base+4, 'MaxEntScore':maxentscore_pot}
                if win_seq[base+1].upper() == 'T':
                    potential_seq = str(win_seq[base-3:base+6])
                    maxentscore_pot = maxent_fast.score5(potential_seq, matrix=matrix5)
                    if maxentscore_pot > 2:
                        splice_dict[potential_seq] = {'Nature':'Donor','Start':var_pos-30+base-3,'End':var_pos-30+base+6, 'MaxEntScore':maxentscore_pot}
    
    if strand == '-':
        for base in range(25,len(win_seq)-26):
            if win_seq[base].upper() == 'G':
                if win_seq[base-1].upper() == 'A':
                    potential_seq = str(win_seq[base-19:base+4])
                    maxentscore_pot = maxent_fast.score3(potential_seq, matrix=matrix3)
                    if maxentscore_pot > 2:
                        splice_dict[potential_seq] = {'Nature':'Acceptor','Start':(len(ref)-var_pos)+win_width-(base-25)+19,'End':(len(ref)-var_pos)+5-(base-25)-3, 'MaxEntScore':maxentscore_pot}
                if win_seq[base+1].upper() == 'T':
                    potential_seq = str(win_seq[base-3:base+6])
                    maxentscore_pot = maxent_fast.score5(potential_seq, matrix=matrix5)
                    if maxentscore_pot > 2:
                        splice_dict[potential_seq] = {'Nature':'Donor','Start':(len(ref)-var_pos)+win_width-(base-25)+3,'End':(len(ref)-var_pos)+5-(base-25)-5, 'MaxEntScore':maxentscore_pot}
    return splice_dict
                

def window_splice_var(var_pos,var_b,strand,ref,win_width):
    if strand == '+':
        return splice_search(ref[var_pos-30:var_pos] + var_b + ref[var_pos+1:var_pos+30+1],var_pos,'+',ref,win_width)
    if strand == '-': 
        var_pos = len(ref)-var_pos
        return splice_search(ref[var_pos-30:var_pos] + var_b + ref[var_pos+1:var_pos+30+1],var_pos,'-',ref,win_width)
      
        
############################################################################
def sub_variant(seq,nature,strand,start,end):
    if nature.upper() == 'DONOR':
        if strand =='-':
            if start-3 == end+5:
                return start-3,seq[3]
        if strand =='+':
            if start+3 == end-5:
                return start+3,seq[3]
    if nature.upper() == 'ACCEPTOR':
        if strand =='-':
            if start-18 == end+4:
                return start-18,seq[18]
        if strand =='+':
            if start+18 == end-5:
                return start+18,seq[18]
        
        
############################################################################
def main(): 
    #Input variant details
    chr_pos = 'chr11'
    var_pos = 47367305
    var_b = 'A'
    var_ref_b = 'G'
    strand = '-'
# =============================================================================
#     chr_pos = 'chr11'
#     var_pos = 47364865
#     var_ref_b = 'C'
#     var_b = 'T'
#     strand = '-'
# =============================================================================
# =============================================================================
#     chr_pos = 'chr1'
#     var_pos = 236906186
#     var_b = 'G'
#     var_ref_b = 'A'
#     strand = '+'
# =============================================================================
# =============================================================================
#     chr_pos = 'chr10'
#     var_pos = 112543080
#     var_b = 'A'
#     var_ref_b = 'G'
#     strand = '+'
# =============================================================================
    
    ref = genome_ref(chr_pos,g,strand)
    
    if strand == '-':
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        var_b = "".join(complement.get(base, base) for base in reversed(var_b))
        var_ref_b = "".join(complement.get(base, base) for base in reversed(var_ref_b))
    
    #check if reference base mach the vcf file
    ref_b = ref[var_pos-1]
    #print(ref_b.upper(),var_ref_b)
    if ref_b.upper() == var_ref_b:
        print('The reference genome nuclotide matches vcf file')
    else:
        print('The reference genome nuclotide cannot match to che reference sequence.')
    
    ##### if a variant change lead to an potential splice site change
    # return the type and sequence
    type_res, seq_res, ref_seq = splice_var(var_pos,var_b,strand,ref)
# =============================================================================
#     print(type_res, seq_res, ref_seq)
#     print(len(seq_res), len(ref_seq))
# =============================================================================
    if type_res.upper() == 'DONOR':
        maxentscore_var = maxent_fast.score5(seq_res, matrix=matrix5)
        if maxentscore_var >= 5:
            maxentscore_ref = maxent_fast.score5(ref_seq, matrix=matrix5)
            maxent_diff = maxentscore_var - maxentscore_ref
            if maxent_diff >= 3:
                print(type_res, maxentscore_var, seq_res)
            else:
                print('This variant is not potential splice donor site.')
        else:
            print('This variant is not potential splice donor site.')
    else:
        maxentscore_var = maxent_fast.score3(seq_res, matrix=matrix3)
        if maxentscore_var >= 5:
            maxentscore_ref = maxent_fast.score3(ref_seq, matrix=matrix3)
            maxent_diff = maxentscore_var - maxentscore_ref
            if maxent_diff >= 3:
                print(type_res, maxentscore_var, seq_res)
            else:
                print('This variant is not potential splice acceptor site.')
        else:
            print('This variant is not potential splice acceptor site.')
    
    #part B
    #####check if there is possible donor/acceptor site in 500bp len 
    width = 500
    adj_var_dict = adj_var_seq(var_pos,var_b,strand,ref,splice_var(var_pos,var_b,strand,ref)[0],width)
    for seq in adj_var_dict:
        #print(len(seq))
        #if str(seq) == 'CTCTACCCTCTTCTGAAAAGAAA':
        #    print('yes')
        if len(seq) <= 15:
            maxentscore_1 = maxent_fast.score5(str(seq), matrix=matrix5)
            if maxentscore_1 >= 2:
                print(seq, maxentscore_1, adj_var_dict[seq])
        else:
            maxentscore_1 = maxent_fast.score3(str(seq), matrix=matrix3)
            if maxentscore_1 >= 2:
                print(seq, maxentscore_1, adj_var_dict[seq])
        
    #part C
    ##### using sliding window method to check if there is possible 'AG' or GT'
    win_width = 5
    win_splice_dict = window_splice_var(var_pos, var_b, strand, ref,win_width)
    for var in win_splice_dict:
        print(var,win_splice_dict[var])
    
    #Calculate sequence score and filter
    #if score > threshold, go to part B
    
    #Part D 
    #for var in win_splice_dict:
    #    sub_var_pos,sub_var_b = sub_variant(var,win_splice_dict[var]['Nature'],strand,win_splice_dict[var]['Start'],win_splice_dict[var]['End'])
    #    sub_adj_var_dict = adj_var_seq(sub_var_pos,sub_var_b,strand,ref,win_splice_dict[var]['Nature'],width)
    #    for seq in sub_adj_var_dict:
    #        print(seq,sub_adj_var_dict[seq])

            
if __name__=="__main__": 
    
    import genomepy
    # package name Biopython
    from Bio.Seq import Seq
    from maxentpy import maxent
    from maxentpy import maxent_fast
    from maxentpy.maxent import load_matrix5, load_matrix3
    
    #Install reference genome
    #genomepy.install_genome('hg19', 'UCSC')
    ## name: hg19
    ## local name: hg19
    ## fasta: /Users/yuchen/.local/share/genomes/hg19/hg19.fa
    ## Created config file /Users/yuchen/Library/Application Support/genomepy/genomepy.yaml
    g = genomepy.Genome('hg19',r"/Users/yuchen/.local/share/genomes")
    
    matrix3 = load_matrix3()
    matrix5 = load_matrix5()

    #main function
    main() 
