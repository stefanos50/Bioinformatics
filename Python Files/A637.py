import itertools
import numpy as np
import pandas as pd
from Bio import SeqIO
m, s, match = [-1,-2,5]

amino_acid_dict = {"A":"Ala","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu","G":"Gly","H":"His","I":"Iie","L":"Leu","K":"Lys"
              ,"M":"Met","F":"Phe","P":"Pro","S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val","R":"Arg"}

#Διαβάζει ένα αρχείο fasta του οποίο το όνομα το δέχεται ως όρισμα
#και επιστρέφει τις ακολουθίες που περιέχει αλλά και τα id
def read_fasta_file(file_name):
    seqs = []
    names = []
    for seq_record1 in SeqIO.parse(file_name, "fasta"):
        names.append(str(seq_record1.id))
        seqs.append(str(seq_record1.seq))
    return seqs,names

#Επιστρέφει την αμινοξικη ακολουθια
def amino_acid_sequence(seq):
    amino_acid = []
    for i in range(len(seq)):
        aseq = ""
        for j in range(len(seq[i])):
            aseq = aseq + amino_acid_dict[seq[i][j]] + '-'
        amino_acid.append(aseq[:-1])
    return amino_acid

#Για μια λίστα w με ακολουθέις επιστρέφει κάθε πιθανό συνδιασμό ανάμεσα
#σε κάθε δυάδα
def w_combinations(w):
    wij = []
    for L in range(0, len(w) + 1):
        for subset in itertools.combinations(w, L):
            if len(subset) == 2:
                wij.append([subset[0],subset[1]])
                wij.append([subset[1],subset[0]])
    return wij

#Τυπώνει τους πίνακες
def print_array(S,T,A):
    col = ' ' + S
    ind = ' ' + T
    df = pd.DataFrame(A.tolist(), columns=[char for char in col], index=[char for char in ind])
    print('\n',df)

#Συνάρτηση που υλοποιεί τη καθολική στοίχιση
def global_alignment(S, T):

    #Αρχικοποίηση του πίνακα με μηδέν
    q = np.full((len(T) + 1, len(S) + 1), ['0'], dtype=str)
    V = np.zeros((len(T) + 1, len(S) + 1))

    #Αρχικοποίηση της γραμμής και στήλης 0 με βάση το gap penalty
    for j in range(1, len(S) + 1):
        V[0][j] = V[0][j-1] + s
        q[0][j] = '←'
    for i in range(1, len(T) + 1):
        V[i][0] = V[i-1][0] + s
        q[i][0] = '↑'

    #Γέμισμα του υπόλοιπου πίνακα
    for i in range(1, len(T) + 1):
        for j in range(1, len(S) + 1):
            horizontal_distance = V[i][j - 1] + s
            vertical_distance = V[i - 1][j] + s
            if T[i - 1] == S[j - 1]:
                diagonal_distance = V[i - 1][j - 1] + match
            else:
                diagonal_distance = V[i - 1][j - 1] + m
            V[i][j] = max(horizontal_distance, vertical_distance, diagonal_distance)
            #Αποθήκευση προέλευσης κάθε τιμής στον πίνακα q
            if max(horizontal_distance, vertical_distance, diagonal_distance) == horizontal_distance:
                q[i][j] = "←"
            elif max(horizontal_distance, vertical_distance, diagonal_distance) == vertical_distance:
                q[i][j] = "↑"
            elif max(horizontal_distance, vertical_distance, diagonal_distance) == diagonal_distance:
                q[i][j] = "↖"
    #Τύπωση των πινάκων V,q
    print_array(S, T, V)
    print_array(S, T, q)

    #Διαδικασία backtracking από τη θέση [-1,-1]
    def backtrack():
        i = len(T)
        j = len(S)
        while i != 0 or j != 0:
            #Αν η τιμή προήλθε διαγώνια τότε τα σύμβολα παραμένουν ως έχουν
            if q[i][j] == "↖":
                i -= 1
                j -= 1
                if T[i] == S[j]:
                    yield T[i], S[j]
                elif T[i] != S[j]:
                    yield T[i], S[j]
            #Αν προήλθε από αριστερά τότε στη πρώτη ακολουθία μπαίνει κενό ενώ η δεύτερη μένει ως έχει
            elif q[i][j] == "←":
                j -= 1
                yield '-', S[j]
            #Αν προήλθε από πάνω τότε στη δεύτερη ακολουθία μπαίνει κενό ενώ η πρώτη μένει ως έχει
            elif q[i][j] == "↑":
                i -= 1
                yield T[i], '-'
            #Αλλιώς κάτι πήγε λάθος
            else:
                assert (False)
    #Επιστροφή της στοίχισης και του score
    return [''.join(reversed(s)) for s in zip(*backtrack())] , V[-1][-1]



sequences,names = read_fasta_file('rcsb_pdb_1BBT.fasta')
amino_acid_sequences = amino_acid_sequence(sequences)
nucleotide_sequences,_ = read_fasta_file('nuc_seq.fasta')

#Εκτύπωση της αμινοξικής και νουκλεοτιδικής ακολουθιας για τα ερωτήματα i,ii
for i in range(len(sequences)):
    print("-------------------------------------------------")
    print("ID ",names[i])
    print("Sequence ",sequences[i])
    print("Amino Acid Sequence ",amino_acid_sequences[i])
    print("Nucleotide Sequence ",nucleotide_sequences[i])
    print("-------------------------------------------------")

global_scores = []
global_allignments = []
v = sequences[0]
sequences.pop(0)
wij = w_combinations(sequences)

#Υπολογισμός καθολικές στοίχισης για κάθε πιθανό συνδιασμό συνένωσης
for i in range(len(wij)):
    concatenation = wij[i][0] + wij[i][1]
    algn,scr = global_alignment(v,concatenation)
    print(algn[0])
    print(algn[1])
    print("With score:",scr)
    global_scores.append(scr)
    global_allignments.append(algn)

#Το index στη λίστα με τα scores που βρίσκεται το max
max_ind = np.argmax(global_scores)

#Τυπώνει τα αποτελέσματα της μέγιστης
print("Chimeric alignment problem:")
print(global_allignments[max_ind][0])
print(global_allignments[max_ind][1])
print("With score: ",global_scores[max_ind])