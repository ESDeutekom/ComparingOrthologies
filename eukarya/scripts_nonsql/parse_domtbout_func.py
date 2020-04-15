#python3
######################################################################
## Parse hmmdomtbout to get non overlapping and best hits from file
## For further processing in parse_domtbout.py
######################################################################
import sys
from operator import itemgetter

def parse_domtbout(domtbfile, accession = False):
    #parse file into a dictionary
    #Pfam name in different columns: difference between accession name and query name
    seq_pfam = {}
    domtbfile = open(domtbfile, "r")
    for line in domtbfile:
        #see which file is being parsed
        if line.startswith("#"):
            if line.rstrip().find('Target file') != -1:
                print(line)
                continue
            continue
        line = line.rstrip().split()
        seq = line[0]
        if accession == True:
            pfam = line[4] #in Julian's data
            pfam = pfam[:pfam.find(".")]
        if accession== False:
            pfam = line[3] #in OG pipeline data
        eval1 = float(line[6]) #sequence e-value
        bit_1 = float(line[7]) # sequence bitscore
        eval2 = float(line[12]) # domain evalue
        bit_2 = float(line[13]) # domain bitscore used to find better hit.
        ali_start = int(line[17]) #alignment start coord
        ali_stop = int(line[18]) #alignment stop coord
        env_start = int(line[19]) #env start coord
        env_stop = int(line[20]) #envstop coord
        if seq in seq_pfam:
            seq_pfam[seq] += [[pfam, eval1, bit_1, eval2, bit_2, ali_start, ali_stop, env_start, env_stop]]
        else:
            seq_pfam[seq] = [[pfam, eval1, bit_1, eval2, bit_2, ali_start, ali_stop, env_start, env_stop]]
    domtbfile.close()
    return seq_pfam

#thought behind this code
#Sort hits per sequences according to best
#then start looking for overlaps from the first best hit
#any overlapping with this best hit will be removed
#then the next best hit (that is not overlapping with the previous best hit) will check for overlap with the next best hits, etc.
def unique_nonoverlaps(hit_dict, overlap = 15, score = "evalue"):
    #takes dictionary with hits
    #finds best non overlapping hits
    for seq, hits in hit_dict.items():
        #sort pfam hits according to highest bit score of domain in descending order
        #In this way, the first one is always the best.
        if score == "bit":
            hits.sort(key = itemgetter(2), reverse=True) #get item for bit score domain
        #Other way around for pvalue. You want the lowest pvalue
        elif score == "evalue":
            hits.sort(key = itemgetter(1), reverse=False)#get item for evalue score and sort lists in list to that value

    #check for overlap in hits in sorted dict
    for (seq, hits) in hit_dict.items():
        i = 0
        print(hits)
        if len(hits) > 1: # only 1 hit means nothing to compare
            hit_l = len(hits[0])
            while i < len(hits): #first to compare
                j = i+1 #always one further from i
                while j < len(hits): #second to compare
                    #pfam, eval1, bit_1, eval2, bit_2, ali_start, ali_stop, env_start, env_stop
                    #itersect gets the numbers that intersect between two lists
                    #The the two lists are the numbers ranging from env_start to env_stop +1
                    intersect = set(list(range(hits[i][-2], hits[i][-1]+1))).\
                    intersection(list(range(hits[j][-2], hits[j][-1]+1))) #check if the range intersects
                    if len(intersect) >= overlap: # the intersect has more # of overlapping parts than overlap
                        hits.pop(j) #remove from list in dictionary
                        j = j #restart counting/stay the same
                    else:
                        j += 1 #go on counting
                i += 1
    for seq, hits in hit_dict.items():
        hits.sort(key = itemgetter(5)) #rearrange hits according to location
    return hit_dict

def get_best_scoring_hits(non_overlap_dict, score = "evalue", score_cutoff = 1e-3):
    #Counts for possible fussions (multiple different parts) and duplications (same part multiple times)
    best_dict = {}
    for seq, hits in non_overlap_dict.items():
        for hit in hits:
            leca_id = hit[0]
            score_now = hit[1]
            if seq not in best_dict:
                best_dict[seq] = {}
                best_dict[seq][leca_id] = [score_now]
            else:
                if leca_id not in best_dict[seq]:
                    if score == "evalue":
                        if score_now <= score_cutoff: #pval smaller to be better
                            best_dict[seq][leca_id] = [score_now]
                    elif score == "bit":
                        if score_now >= score_cutoff: #pval smaller to be better
                            best_dict[seq][leca_id] = [score_now]
                else:
                    if score == "evalue":
                        if score_now <= score_cutoff:
                            best_dict[seq][leca_id] += [score_now]
                    elif score == "bit":
                        if score_now >= score_cutoff:
                            best_dict[seq][leca_id] += [score_now]
    return best_dict

def transpose_hits(best_hits):
    transposed_hits = {}
    for seq, hit_list in best_hits.items():
        for hit in hit_list:
            if hit in transposed_hits:
                transposed_hits[hit] += [seq]
            else:
                transposed_hits[hit] = [seq]
    return transposed_hits
