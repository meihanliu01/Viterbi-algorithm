import pyfastx
import gzip
import pandas as pd
import numpy as np
import math
from collections import defaultdict
import gffpandas.gffpandas as gffpd
####### a
# files
file1 = 'Vibrio_cholerae.GFC_11.dna.toplevel.fa'
file2 = 'Vibrio_cholerae.GFC_11.37.gff3.gz'
file3 = 'Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa'
file4 = 'Vibrio_vulnificus.ASM74310v1.37.gff3'
#read fa file
def readFa(file):
    names = []
    seqs = []
    for name, seq in pyfastx.Fasta(file, build_index=False):
        names.append(name)
        seqs.append(seq)
    return names, seqs
#read gff3 file
in_handle = gzip.open(file2, "rt")
lines = in_handle.readlines()
annotated_gene_list = []
annotated_gene_list.append([])
contig = 'DN38.contig00001'
num = 0 # number of contig sequence
for line in lines:
    if ("supercontig" in line):
        if (contig not in line):
            contig = line.split()[0]
            num += 1
            annotated_gene_list.append([])
    if (("CDS" in line) and ("+" in line)):
        line_array = line.split()
        start = int(line_array[3])
        end = int(line_array[4])
        annotated_gene_list[num].append([start, end])
in_handle.close()
# read fa file: file1
names,seqs = readFa(file1)
seqs[89],seqs[88] = seqs[88], seqs[89] # the order of sequence in file1 is wrong
seqs[116],seqs[117] = seqs[117], seqs[116]
total_nucleotides_num = 0
for i in range(0, len(seqs)):
    total_nucleotides_num = total_nucleotides_num + len(seqs[i])
#print("total_nucleotides_num: " + str(total_nucleotides_num))

#genic_len: number of nt in all genic regions
#genic_num: number of genic regions
genic_len, genic_num = 0,0
intergenic_len, intergenic_num = 0,0
genic_list = [] #string list
intergenic_list = []
last_end = -1
genic_A, genic_C, genic_T, genic_G = 0,0,0,0
total_A, total_C, total_T, total_G = 0,0,0,0

for i in range(0, len(seqs)):
    if annotated_gene_list[i] != []:
        for start, end in annotated_gene_list[i]:
            genic_list.append(seqs[i][start : end + 1])
            genic_len += end - start + 1
            for k in seqs[i][start:end + 1]:
                if k == 'A':
                    genic_A += 1
                if k == 'T':
                    genic_T += 1
                if k == 'C':
                    genic_C += 1
                if k == 'G':
                    genic_G += 1
    for j in seqs[i]:
        if i == 'A':
            total_A += 1
        if i == 'T':
            total_T += 1
        if i == 'C':
            total_C += 1
        if i == 'G':
            total_G += 1
genic_avg_len = genic_len/len(genic_list)
intergenic_avg_len = (total_nucleotides_num-genic_len)/(len(genic_list)+1)
# print("genic_avg_len: "+ str(genic_avg_len))
# print("intergenic_avg_len: " + str(intergenic_avg_len))
intergenic_A = total_A - genic_A
intergenic_T = total_T - genic_T
intergenic_C = total_C - genic_C
intergenic_G = total_G - genic_G
intergenic_len = intergenic_A + intergenic_C + intergenic_G + intergenic_T
freq_A = intergenic_A/intergenic_len
freq_C = intergenic_C/intergenic_len
freq_T = intergenic_T/intergenic_len
freq_G = intergenic_G/intergenic_len
lst1 = ["A", "T", "C", "G"]
lst2 = [freq_A, freq_C, freq_G, freq_T]
nucleotide_freq_tbl = list(zip(lst1, lst2))
#print("nucleotide_freq_tbl: " + str(nucleotide_freq_tbl))
nucleotide_freq_dict = {'A': freq_A, 'C': freq_C, 'T': freq_T, 'G': freq_G}
# print(nucleotide_freq_dict)

dict_codons={}
codon_num = 0
dict_start={}
dict_end={}
start_num, end_num = 0,0

for i in range(0, len(seqs)):
    if annotated_gene_list[i] != []:
        for start, end in annotated_gene_list[i]:
            start_seq = seqs[i][start-1:start+2]
            end_seq = seqs[i][end-3:end]
            if start_seq in dict_start:
                dict_start[start_seq]+=1
                start_num +=1
            else:
                dict_start[start_seq] = 1
                start_num += 1
            if end_seq in dict_end:
                dict_end[end_seq]+=1
                end_num+=1
            else:
                dict_end[end_seq] = 1
                end_num += 1
            interval = (seqs[i][n:n + 3] for n in range(start + 2, end - 4, 3))
            for codon in interval:
                if codon in dict_codons:
                    dict_codons[codon] += 1
                else:
                    dict_codons[codon] = 1
            codon_num += int((end - start + 1)/3)

for items, values in dict_start.items():
    dict_start[items] = values/start_num
for items, values in dict_end.items():
    dict_end[items] = values/end_num
for items, values in dict_codons.items():
    dict_codons[items] = values/codon_num

dict_codons.update([('TAA',0),('TGA',0),('TAG',0)])
#print(dict_codons)

####### part b
#######Write a Fasta file
#######

####data may be useful
avg_genic_len = 985
ave_intergenic_len = 1134
# number of nt in start codons
start_nt = 3 * start_num
end_nt = 3 * end_num

###### wirte emission matrix
######
###### initial state probability
initial_prob_dict = {'Intergenic': 1, 'Start': 0, 'Middle': 0, 'Stop': 0}
#print(initial_prob)
###### transition probabilities matrix
transition_dict = {'Intergenic': {'Intergenic': 1133/1134, 'Start': 1/1134, 'Middle': 0, 'Stop': 0},
                   'Start': {'Intergenic': 0, 'Start': 0, 'Middle': 1, 'Stop': 0},
                   'Middle': {'Intergenic': 0, 'Start': 0, 'Middle': 985/988, 'Stop': 3/988},
                   'Stop': {'Intergenic': 1, 'Start': 0, 'Middle': 0, 'Stop': 0}
                   }
###### emission state probability
inter_emission = {}
start_emission = {}
mid_emission = dict_codons
stop_emission = {}
for key, value in dict_codons.items():
    if key == 'ATG':
        start_emission[key] = dict_start.get(key)
    elif key == 'GTG':
        start_emission[key] = dict_start.get(key)
    elif key == 'TTG':
        start_emission[key] = dict_start.get(key)
    else:
        start_emission[key] = 0
#print(start_emission)
# start_emission = dict_start
for key, value in dict_codons.items():
    if key == 'TAA':
        stop_emission[key] = dict_end.get(key)
    elif key == 'TGA':
        stop_emission[key] = dict_end.get(key)
    elif key == 'TAG':
        stop_emission[key] = dict_end.get(key)
    else:
        stop_emission[key] = 0
inter_emission = nucleotide_freq_dict
#write config file
f = open("input.txt", "w")
f.write("Transtion matrix: " + '\n')
f.write(str(transition_dict) + '\n')
f.write("Intergenic emission: " + '\n')
f.write(str(inter_emission) + '\n')
f.write("Start emission: " + '\n')
f.write(str(start_emission) + '\n')
f.write("Middle emission: " + '\n')
f.write(str(mid_emission)+ '\n')
f.write("Stop emission: " + '\n')
f.write(str(stop_emission)+ '\n')
f.write("Initial probability table: " + '\n')
f.write(str(initial_prob_dict) + '\n')
f.write("Average length of intergenic regions: " + '\n')
f.write(str(intergenic_avg_len) + '\n')
f.write("Average length of genic region: " + '\n')
f.write(str(genic_avg_len) + '\n')
f.write("The nucleotide frequency table for intergenic regions: " + '\n')
f.write(str(nucleotide_freq_tbl) + '\n')
f.write("The codon frequency table for genic regions: " + '\n')
f.write(str(dict_codons) + '\n')
f.close()
# viterbi algo
def viterbi(file = file3, infilename = 'input.txt'):
    f = open(infilename, 'r').read().splitlines()
    transition_dict, inter_emit, start_emit, mid_emit,stop_emit, initial_prob  = eval(f[1]), eval(f[3]), eval(f[5]), eval(f[7]), eval(f[9]), eval(f[11])
    # convert all probabilities to negative log
    def convertLog(dict):
        new_log_matrix = {}
        for key, value in dict.items():
            if value != 0:
                new_log_matrix[key] = -np.log(value)
            else:
                new_log_matrix[key] = math.inf
        return new_log_matrix
    log_transition_matrix=transition_dict
    for k in transition_dict.keys():
        for key, value in transition_dict[k].items():
            if transition_dict[k][key] != 0:
                log_transition_matrix[k][key] = -np.log(value)
            else:
                log_transition_matrix[k][key] = math.inf
    log_initial_matrix = convertLog(initial_prob)
    log_inter_emission = convertLog(inter_emit)
    log_start_emission = convertLog(start_emit)
    log_stop_emission = convertLog(stop_emit)
    log_mid_emission = convertLog(mid_emit)
    # read input file
    names, seqs = readFa(file)
    file_ = pyfastx.Fasta(file)
    uni_names = np.unique(names)
    out = open('output.gff3', 'w') # output file
    for gen in uni_names:
        seq = file_[str(gen)]
        prob = np.full((4, len(seq)), math.inf) # a matrix store probabilities
        check = np.full((4, len(seq)), False) # check whether the current nt i can be the first nt of a codon (i:i+3)
        prev = np.full((4, len(seq)), math.inf) # a matrix store the previous states
        prob[0,0] = log_initial_matrix['Intergenic'] + log_inter_emission[seq[0]] #initialize
        for i in range(1,len(seq)):
            if log_inter_emission[seq[i]] == math.inf: #not inter
                prob[0,i] = math.inf
                break
            I = prob[0,i-1] + log_inter_emission[seq[i]] + log_transition_matrix['Intergenic']['Intergenic']
            SP = math.inf
            if i-3 >= 0 and check[3,i-3]:
                SP = prob[3,i-3] + log_inter_emission[seq[i]] + log_transition_matrix['Stop']['Intergenic']
            prob[0, i] = min(I, SP)
            if prob[0,i] != math.inf: #store the prev state
                if prob[0,i] == I:
                    prev[0, i] = 0
                else:
                    prev[0, i] = 3
            if i+3 < len(seq):
                c = str(seq[i:i+3]) #codon
                if log_start_emission[c] == math.inf: # not start
                    prob[1,i] = math.inf
                else:
                    I = prob[0,i-1] + log_start_emission[c] + log_transition_matrix['Intergenic']['Start']
                    prob[1,i:i+3] = I # do i need to check if p0 is inf??????
                    prev[1,i] = 0 # intergenic->start
                    prev[1,i+1], prev[1,i+2], prev[1,i+3] = 1,1,1 # start codon
                    check[1, i] = True
                if log_stop_emission[c] == math.inf: #not stop
                    prob[3,i] = math.inf
                else:
                    M = math.inf
                    if i-3 >= 0 and check[2,i-3]:
                        M = prob[2,i-3] + log_stop_emission[c] + log_transition_matrix['Middle']['Stop'] # only stop -> middle can happen
                    if M != math.inf:
                        prob[3,i:i+3] = M
                        check[3,i] = True
                        prev[3,i] = 2
                        prev[1,i+1], prev[1,i+2], prev[1,i+3] = 1,1,1
                    else:
                        prob[3,i] = math.inf
                if log_mid_emission[c] == math.inf: #not middle
                    prob[2,i] = math.inf
                else:
                    S, M = math.inf, math.inf
                    if i-3 >= 0 and check[1,i-3]:
                        S = prob[1,i-3] + log_mid_emission[c] + log_transition_matrix['Start']['Middle']
                    if i-3 >= 0 and check[2,i-3]:
                        M = prob[2,i-3] + log_mid_emission[c] + log_transition_matrix['Middle']['Middle']
                    m = min(S, M)
                    if m != math.inf:
                        prob[2,i:i+3] = m
                        check[2,i] = True
                        if m == S:
                            prev[2,i] = 1
                        if m == M:
                            prev[2,i] = 2
                        prev[2, i + 1], prev[2, i + 2], prev[2, i + 3] = 2,2,2
                    else:
                        prob[2,i] = math.inf
        #traceback
        state = np.empty(len(seq), dtype=int) # backtrack the states for every nt
        state[len(seq)-1] = np.argmin(prob[:,len(seq)-1])
        j = len(seq)-2
        k = 0
        while j >= 0:
            state[j] = prev[int(state[j+1]), j+1]
            if state[j] == 1 or state[j] == 2 or state[j] == 3:
                state[j-2:j] = state[j]
                j -= 3
            else:
                j -= 1
        while k < len(state):
            if state[k] == 1:
                start = k+1
                while state[k] != 0:
                    k += 1
                end = k
                out.write(f'{gen}\tena\tCDS\t{start}\t{end}\t.\t+\t0\n')
            k += 1
    out.close()
viterbi(file3, 'input.txt')


########## d and e
def match(file1, file2): # file1: file4
    names, seqs = readFa('Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa')
    f1 = gffpd.read_gff3(file1)
    f1_ = f1.df
    filter1 = f1_.loc[(f1_['strand'] == '+') & (f1_['type'] == 'CDS')]
    id_file1 = np.unique(filter1['seq_id'])
    f2 = gffpd.read_gff3(file2)
    f2_ = f2.df
    filter2 = f2_.loc[(f2_['strand'] == '+') & (f2_['type'] == 'CDS')]
    both, start, end = 0,0,0
    total_num = len(filter1)
    total_num2 = len(filter2)
    perfect = dict()
    for id in id_file1:
        seq_id1 = filter1.loc[filter1['seq_id'] == id]
        seq_id2 = filter2.loc[filter2['seq_id'] == id]
        for i in range(0, seq_id1.shape[0]):
            for j in range(0, seq_id2.shape[0]):
                if seq_id1.iloc[i]['start'] == seq_id2.iloc[j]['start'] and seq_id1.iloc[i]['end'] == seq_id2.iloc[j]['end']:
                    both += 1
                    if id not in perfect:
                        perfect[id] = []
                    perfect[id].append([int(seq_id1.iloc[i]['start']), int(seq_id1.iloc[i]['end'])])
                if seq_id1.iloc[i]['start'] == seq_id2.iloc[j]['start'] and seq_id1.iloc[i]['end'] != seq_id2.iloc[j]['end']:
                    start += 1
                if seq_id1.iloc[i]['start'] != seq_id2.iloc[j]['start'] and seq_id1.iloc[i]['end'] == seq_id2.iloc[j]['end']:
                    end += 1
    both_ratio = both / total_num
    start_ratio = start / total_num
    end_ratio = end / total_num
    neither_ratio = (total_num - both - start - end) / total_num

    both_ratio2 = both / total_num2
    start_ratio2 = start / total_num2
    end_ratio2 = end / total_num2
    neither_ratio2 = (total_num2 - both - start - end) / total_num2

    print('The fraction of annotated genes on the positive strand:')
    print('both: '+ str(both_ratio))
    print('start: ' + str(start_ratio))
    print('end: ' + str(end_ratio))
    print('neither: ' + str(neither_ratio))
    print('The fraction of your predicted genes:')
    print('both: ' + str(both_ratio2))
    print('start: ' + str(start_ratio2))
    print('end: ' + str(end_ratio2))
    print('neither: ' + str(neither_ratio2))

    len_perfect = 0
    sum_perfect = 0
    for key,value in perfect.items():
        len_perfect += len(perfect[key])
        for i in value:
            sum_perfect += i[1] - i[0] +1
    print(sum_perfect/len_perfect)

    def findIndex(input):
        index = 0
        for name in names:
            if name == input:
                return index
            index += 1
    start_codon = []
    stop_codon = []
    for contig,list in perfect.items():
        contig_input = findIndex(contig)
        for pair in list:
            start = pair[0]
            end = pair[1]
            start_codon.append(seqs[contig_input][start:start + 3])
            stop_codon.append(seqs[contig_input][end - 3:end])
    s,e={},{}
    for i in start_codon:
        if i in s:
            s[i] += 1
        else:
            s[i] = 1
    print(s)
    for j in stop_codon:
        if j in e:
            e[j] += 1
        else:
            e[j] = 1
    print(e)

match(file4, 'output.gff3')


