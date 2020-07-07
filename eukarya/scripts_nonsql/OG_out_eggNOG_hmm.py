 #!/hosts/linuxhome/scarab/eva2/Programs/miniconda3/bin/python
#python3
######################################################
## parse eggNOG hmmer file to retrieve best OG hits
######################################################
import sys
import gzip

annotation_in=sys.argv[1]
ortho_out = sys.argv[2]

annotation_in = gzip.open(annotation_in,'rt')

hit_dict = {}
#eggNOG mapper with best hit annotation for euNOG in 11th column
for line in annotation_in:
    line.rstrip()  # Remove newline and trailing spaces
    print(line)
    if line[0] == "#":
        continue
    elements = line.rsplit("\t")
    protein_id = elements[0] #eukarya4 name
    best_og_model = elements[10].split("|") #best_eugNOG|evalue|score
    OG = best_og_model[0]
    e_value = best_og_model[1]
    euNOG = "".join([OG,'@euNOG'])
    #Anything passing an e-value of 10-3 is considered a hit
    #seem eggNOG mapper did this already, but just in case
    if float(e_value) <= 1e-3:
        if euNOG in hit_dict:
            hit_dict[euNOG] += [protein_id]
        else:
            hit_dict[euNOG] = [protein_id]
annotation_in.close()
ortho_out = open(ortho_out, 'w')
for key, value in hit_dict.items():
    ortho_out.write(" ".join([key+':', " ".join(value), "\n"]))
