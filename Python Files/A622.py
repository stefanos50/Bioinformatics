import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
#m = mismatch penalty
#s = gap penalty
#match = match reward
m, s, match = [-1,-5,4]

#Επιστρέφει το index του κελιού που έχει
#τη μεγαλύτερη τιμή στη τελευταία γραμμή
def get_max(row):
    row.pop(0)
    return np.argmax(row) + 1

#Τυπώνει τους πίνακες
def print_array(S,T,A):
    col = ' ' + S
    ind = ' ' + T
    df = pd.DataFrame(A.tolist(), columns=[char for char in col], index=[char for char in ind])
    print('\n',df)

#Διαβάζει fasta αρχείο και επιστέφει την ακολουθία του
def read_fasta_file(file_name):
    for seq_record1 in SeqIO.parse(file_name, "fasta"):
        sequence = seq_record1.seq
    return sequence

def overlap_alignment(S,T):
        q = np.full((len(T)+1,len(S)+1), ['0'],dtype=str)
        V = np.zeros((len(T)+1,len(S)+1))
        #Από το κελί 1,1 και για κάθε άλλο κελί
        for i in range(1,len(T)+1):
            for j in range(1,len(S)+1):
                #Περίπτωση προσθαφαιρέσεων
                horizontal_distance = V[i][j-1] + s
                vertical_distance = V[i-1][j] + s
                #Αν είναι το ίδιο σύμβολο και στις δύο ακολουθίες
                #τότε επιβραβεύεται
                if T[i-1] == S[j-1]:
                    diagonal_distance = V[i-1][j-1] + match
                #Αν όχι προστήθεται το penalty για ασυμφωνία
                else:
                    diagonal_distance = V[i-1][j-1] + m
                #Το μεγαλύτερο από τα παραπάνω θα μπει στο κελί i,j
                V[i][j] = max(horizontal_distance,vertical_distance,diagonal_distance)
                #Στον πίνακα q αποθηκεύεται από που προήλθε η κάθε τιμή
                if max(horizontal_distance,vertical_distance,diagonal_distance) == horizontal_distance:
                    q[i][j] = "←"
                elif max(horizontal_distance,vertical_distance,diagonal_distance) == vertical_distance:
                    q[i][j] = "↑"
                elif max(horizontal_distance,vertical_distance,diagonal_distance) == diagonal_distance:
                    q[i][j] = "↖"
        print_array(S,T,V)
        print_array(S,T,q)
        def backtrack():
                i = len(T)
                j = len(S)
                max_j = get_max(V[i].tolist())
                reach_max = False
                while i != 0 or j != 0:
                    if j == max_j:
                        reach_max = True
                    #Αν πάμε διαγώνια τότε μένουν όπως έχουν
                    if q[i][j] == "↖" and reach_max == True:
                        i -= 1
                        j -= 1
                        if T[i] == S[j]:
                            yield T[i], S[j]
                        elif T[i] != S[j]:
                            yield T[i], S[j]
                    #Αν πάμε οριζόντια τότε μπαίνει - στην πρώτη ακολουθία
                    elif q[i][j] == "←" and reach_max == True:
                        j -= 1
                        yield '-', S[j]
                    #Αν πάμε πάνω τότε μπαίνει - στην δεύτερη ακολουθία
                    elif q[i][j] == "↑" and reach_max == True:
                        i -= 1
                        yield T[i], '-'
                    #Στην τελευταία γραμμή μέχρι να φτάσουμε από τη θέση -1,-1 στο max
                    #μπαίνει - στη πρώτη ακολουθία
                    elif reach_max == False:
                        j = j - 1
                        yield '-', S[j]
                    #Αν φτάσει στις στήλες i=0 ή j=0 οι οποίες είναι μηδενικές
                    #τότε πηγαίνει προς τα πάνω ή προς τα αριστερά προσθέτοντας
                    # - στην πρώτη ακολουθία ή στη δεύτερη αντίστοιχα
                    elif i==0:
                        j = j - 1
                        yield '-', S[j]
                    elif j==0:
                        i = i -1
                        yield T[i], '-'
                    #Αλλιώς κάτι έχει πάει λάθος
                    else:
                       assert(False)
        #Επιστρέφει πίνακα διαστάσεων 1x2 με την στοίχιση των δύο ακολουθιών
        return [''.join(reversed(s)) for s in zip(*backtrack())]



seq1 = read_fasta_file('αlactalb.fasta')
seq2 = read_fasta_file('c-typelysozyme.fasta')
result = overlap_alignment(seq1,seq2)
print("\n")
print(result[1])
print(result[0])