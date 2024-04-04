import gzip
import matplotlib.pyplot as plt
from pprint import pprint
import numpy as np
import math
from sklearn.decomposition import PCA
from scipy.stats import chi2_contingency
from itertools import combinations

def gc_counter(seq):
    return seq.count('C') + seq.count('G')

def base_counter(seq):
    return seq.count('G') + seq.count('C') + seq.count('A') + seq.count('T')


def gc_content_counter(seq):
    g_counter = seq.count('G')
    c_counter = seq.count('C')
    a_counter = seq.count('A')
    t_counter = seq.count('T')
    return (g_counter + c_counter) / (g_counter + c_counter + a_counter + t_counter)

def parsing_window_1(seq,start,end,pace):
    window = []
    for i in range(len(seq)):
        while end <= i: 
            window.append(seq[start:end:pace])
            start += pace
            end += pace
    return window


def parsing_window_2(seq,start,end,pace):
    window = []
    for i in range(len(seq)):
        while end <= i:
            window.append(seq[start:end:pace])
            start += 1
            end += 1
    return window


def seq_generator(name):
    with gzip.open(name, 'rt') as f:
        dic = {}
        line = f.readline()
        in_first_line = True
        while line != '':
            if line[0] == '>':
                if in_first_line:
                    in_first_line = False
                    dic["header"] = line[:-1]
                    seq = ""
                else:
                    dic["sequence"] = seq
                    yield dic
 
                    dic = {}
                    dic["header"] = line[:-1]
                    seq = ""
            else:
                seq += line[:-1]
 
            line = f.readline()
        else:
            if in_first_line == True:
                yield None
            else:
                dic["sequence"] = seq
                yield dic

def get_I(seq, k):
    
    xy_list = []
    gc_counts_list = []
    
    for pairs in parsing_window_2(seq, 0, k + 1, k):
        gc_counts_list.append(gc_counter(pairs))
        xy_list.append(pairs)

    gc_counts = 0 
    for i in gc_counts_list:
        gc_counts += i 
        

    Nk = {'AA': 0, 'AT': 0, 'AC': 0, 'AG': 0, 'CC': 0, 'CT': 0, 'CA': 0,'CG': 0, 'GG': 0,'GT': 0, 'GA': 0, 'GC': 0,'TT': 0,'TA': 0, 'TC':0, 'TG':0}
    for pair in xy_list:
        Nk[pair] += 1

    Sk = 0
    for count in Nk.values():
        Sk += count

    Pk = {'AA': 0, 'AT': 0, 'AC': 0, 'AG': 0, 'CC': 0, 'CT': 0, 'CA': 0,'CG': 0, 'GG': 0,'GT': 0, 'GA': 0, 'GC': 0,'TT': 0,'TA': 0, 'TC':0, 'TG':0}
    for pair in Pk:
        Pk[pair] = Nk[pair] / Sk

    P = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base in seq:
        P[base] += 1 

    seq_len = len(seq)
    for base in P:
        P[base] = P[base] / seq_len


    Ik = 0
    for x in ['A', 'T', 'G', 'C']:
        for y in ['A', 'T', 'G', 'C']:
            Ik += Pk[x + y] * (math.log10(Pk[x + y] / (P[x] * P[y])))

    return Ik

# compute the chi-squared distance using above formula 
def chi2_distance(a, b): 
  
    chi2, p, dof, ex = chi2_contingency(np.vstack((a, b)))

    #chi =  np.sum([((a - b) ** 2) / (a + b) for (a, b) in zip(a, b)])
    
    return chi2


def part_1():

    seq_1 = next(seq_generator('covid_fasta.gz'))['sequence']

    gc_content_list = []
    window_starts = []
    window_start = 0

    print(len(seq_1))
    print(len(parsing_window_1(seq_1, 0, 500, 30)))
    for pair in parsing_window_1(seq_1, 0, 500, 30): 
        gc_content_list.append(gc_content_counter(pair))
        window_starts.append(window_start)
        window_start += 30 

    x = window_starts
    y = gc_content_list
    plt.plot(x,y, c = 'magenta', linewidth = 2)

    plt.show()


def part_2():
    
    gen = seq_generator('covid_fasta.gz')
    next(gen)
    
    second_strain = next(gen)
    seq = second_strain['sequence']

    ik_list = []
    for k in range(5,455):
        ik_list.append(get_I(seq, k))


    x = [i for i in range(5,455)] 
    y = ik_list
   

    plt.plot(x,y, c = 'magenta', linewidth = 2)

    plt.show()


def part_3():

    gen = seq_generator('covid_fasta.gz')
    first_iteration = True
    i = 0 

    while i < 10:
        
        try:
            seq = next(gen)['sequence']
        except StopIteration:
            print("STAMATHSE")
            break

        k = len(seq) - 27000
        p1_list = parsing_window_1(seq, 0, k, 30)
    
        gc_content_list = []
        for subseq in p1_list:
            gc_content_list.append(gc_content_counter(subseq))
        gc_nparray = np.array(gc_content_list)
        
        ik_list = []
        for k in range(5,455):
            ik_list.append(get_I(seq, k))
        ik_nparray = np.array(ik_list)

        if first_iteration:
            result_array1 = gc_nparray
            result_array2 = ik_nparray
            first_iteration = False
        else:
            result_array1 = np.vstack((result_array1, gc_nparray))
            result_array2 = np.vstack((result_array2, ik_nparray))
        i += 1
        
    return (result_array1,result_array2)


def part_4():
    total_arrays = part_3()
    
    pca1 = PCA(n_components=2)
    pca2 = PCA(n_components = 2)
    pca1.fit((total_arrays)[0])
    pca2.fit((total_arrays)[1])

    array1_PCA = pca1.transform((total_arrays)[0])
    array2_PCA = pca2.transform((total_arrays)[1])

    fig, ax1 = plt.subplots()
    fig, ax2 = plt.subplots()
    
    ax1.plot(array1_PCA[:100,0], array1_PCA[:100,1],'b.')
    ax2.plot(array2_PCA[:100,0], array2_PCA[:100,1],'r.')

    plt.show()
    
    
def part_5():
    gen = seq_generator('covid_fasta.gz')
    count_1 = 0
    first_iteration = True
    count = 0
    
    while True:
        print("compiling #%s" % count)
        try:
            seq = next(gen)['sequence']
        except StopIteration:
            break
        
        codons = {'TTT': 0, 'TGT': 0, 'TAA': 0, 'TGA': 0, 'AAA': 0, 'TAT': 0, 'AGA': 0, 'TGC': 0, 'TGG': 0, 'AAT': 0, 'ACA': 0, 'GTT': 0, 'TTA': 0, 'CAA': 0, 'TAC': 0, 'ACT': 0, 'AAC': 0, 'TTC': 0, 'ATT': 0, 'CTT': 0, 'AGT': 0, 'TCT': 0, 'GAA': 0, 'TTG': 0, 'TCA': 0, 'CAT': 0, 'GCT': 0, 'GAT': 0, 'CAC': 0, 'ACC': 0, 'GGT': 0, 'AGC': 0, 'TAG': 0, 'AGG': 0, 'GAC': 0, 'GTA': 0, 'ATG': 0, 'GCA': 0, 'ATA': 0, 'CCA': 0, 'ATC': 0, 'GGA': 0, 'AAG': 0, 'CTA': 0, 'CCT': 0, 'GTC': 0, 'GTG': 0, 'CAG': 0, 'GAG': 0, 'GGC': 0, 'CTC': 0, 'CTG': 0, 'TCC': 0, 'GCC': 0, 'ACG': 0, 'CGT': 0, 'GGG': 0, 'CCC': 0, 'CGC': 0, 'TCG': 0, 'CGG': 0, 'CGA': 0, 'CCG': 0, 'GCG': 0}

        codons_list = []
        for i in range(0,len(seq),3):
            codon = seq[i:i+3]   
            codons_list.append(codon)     

        for codon in codons_list:
            if len(codon) != 3:
                pass
            else:
                codons[codon] += 1
        
        codons_nparray = np.array(list(codons.values()))

        if first_iteration:
            total_codons_array = codons_nparray
            first_iteration = False       
        else:
            total_codons_array = np.vstack((total_codons_array,codons_nparray))

        count += 1

    to_dump = []
    result = []
    count = 0
    for pivot_list in total_codons_array:
        print("calculating #%s" % count)
        for other_list in total_codons_array:
            
            if chi2_distance(pivot_list, other_list) < 1:
                to_dump.append(count)
                break

        count += 1

    count = 0
    for l in total_codons_array:
        if count not in to_dump:
            result.append(l)
        count += 1

    print(len(result))
            

part_5()