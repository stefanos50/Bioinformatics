import numpy as np
from Bio import SeqIO
import pandas as pd

#Τυπώνει ποιος παίκτης νίκηση ανάλογα το αποτέλεσμα της
#θέσης -1,-1 του πίνακα στρατιγικής
def winning_player(n,m,res):
    if res == 'W':
        print("For sequences with n =",len(n),"and m =",len(m),"the player who will win is player_1")
    else:
        print("For sequences with n=",len(n),"and m=",len(m),"the player who will win is player_2")

#Διαβάζει fasta αρχείο και επιστέφει την ακολουθία του
def read_fasta_file(file_name):
    for seq_record1 in SeqIO.parse(file_name, "fasta"):
        sequence = seq_record1.seq
    return sequence

#Εισάγει στον πίνακα στρατιγικής της θέσης οι οποίες
#είναι γνωστές ότι νικάει ο παίκτης 1
def known_winning_pos(M,n,m):
    for i in range(m+1):
        M[i][0] = "W"
    for j in range(n+1):
        M[0][j] = "W"
    if m>0 and n>0:
        M[1][1] = "W"
    if m >= 2 and n != 0:
        for i in range(2,m+1):
            M[i][1] = "L"
    if n>=2 and m != 0:
        for j in range(2,n+1):
            M[1][j] = "L"
    return M

#Τυπώνει τον πίνακα στρατιγικής όπως στο βιβλίο
def print_array(n,m,A):
    df = pd.DataFrame(A.tolist(), columns=[j for j in range(len(n)+1)], index=[i for i in range(len(m)+1)])
    print('\n',df,'\n')

#Δεδομένου μιας ακολουθίας μήκους n
#και μιας δεύτερης μήκους m
#υπολογίζει τον πίνακα στρατηγικής
#και επιστρέφει τη θέση -1,-1 του πίνακα
#δηλαδή W αν κερδίζει ο παίκτης 1 ή L αν χάνει
def strategy_matrice(n,m):
    R = []
    for i in range(len(m) + 1):
        R.append([0] * (len(n) + 1))
    R = known_winning_pos(R,len(n),len(m))
    if len(n) >= 2 and len(m) >=2:
        for i in range(2,len(m)+1):
            for j in range(2,len(n)+1):
                if R[i-2][j-1] == "W" and R[i-1][j-2] == "W":
                    R[i][j] = 'L'
                else:
                    R[i][j] = 'W'
    print("The strategy matrix is:")
    print_array(n,m,np.array(R))
    return R[-1][-1]


seq1 = read_fasta_file('liver.fasta')
seq2 = read_fasta_file('muscle.fasta')
result = strategy_matrice(seq1,seq2)
winning_player(seq1,seq2,result)