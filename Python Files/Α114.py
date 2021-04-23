import sys
import numpy as np
import pandas as pd

#Q = ονόματα καταστάσεων
Q = ['a','b']
#A = Ονόματα συμβόλων αλληλουχιών
A = ['A','G','T','C']
#N = Αριθμιτική αναπαράσταση συμβόλων
N = [1,2,3,4]

#Τυπώνει τους πίνακες σχέσεων του HMM
def print_HMM(p,t,e):
    df = pd.DataFrame(p.tolist(), columns=['Start'], index=['State a', 'State b'])
    print(df)
    df = pd.DataFrame(t.tolist(), columns=['State a', 'State b'], index=['State a', 'State b'])
    print(df)
    df = pd.DataFrame(e.tolist(), columns=['A', 'G', 'T', 'C'], index=['State a', 'State b'])
    print(df)

#Μετατρέπει την αλληλουχία σε αριθμιτική λίστα όπου
#Α = 1 , G = 2 , T = 3 , C = 4
def sequence_str_to_numeric_list(sequence):
    list = []
    for i in range(len(sequence)):
        list.append(N[A.index(sequence[i])])
    return list

#Με όρισμα μια λίστα επιστρέφει πόσες φορές εμφανίζεται το max element της λίστας
def count_max(prob):
    max_ = max(prob)
    times = prob.count(max_)
    if times > 1:
        return True
    else:
        return False

#Στην περίπτωση που υπάρχουν πολλές ισοπίθανες καταστάσεις με τη max πιθανότητα
#επιστρέφει το index της κατάστασης που θα επιλεγεί σε αυτή τη περίπτωση
def multiple_max_possibility(prob,trans,kat,observation):
    max_ = max(prob)
    index = -sys.maxsize - 1
    p = float('-inf')
    for rows in range(len(prob)):
        if prob[rows] == max_:
            if trans[kat][rows] > p:
                p = trans[kat][rows]
                index = rows
    return index

#Μετατρέπει τη τελική ακολουθία καταστάσεων από λίστα σε string
def states_list_to_str(states_path):
    str = ""
    for char in range(len(states_path)):
        if char == len(states_path)-1:
            str += states_path[char]
            break
        str = str + states_path[char] + " → "
    return str

#Υλοποιεί τον αλγόριθμο viterbi
def viderbi_algorithm(pi,trans,emis,sequence):
    num_hid = trans.shape[0]  # number of hidden states
    num_obs = len(sequence)  # number of observations (not observation *states*)
    print("Number of hidden states",num_hid)
    print("Number of observations",num_obs)
    V = np.zeros((num_hid, num_obs))
    for i in range(num_obs):
        for j in range(num_hid):
            if i == 0:
                V[j][i] = pi[j] + emis[j][sequence[i]-1]
            else:
                p = []
                for k in range(num_hid):
                    p.append(V[k][i-1]+trans[k][j])
                V[j][i] = emis[j][sequence[i]-1] + max(p)
    V = np.array(V)
    df = pd.DataFrame(V.tolist(), columns=[A[char-1] for char in sequence], index=Q)
    print("\n",df,"\n")
    katastaseis = []
    for column in range(num_obs):
        probabilities = []
        for row in range(num_hid):
            probabilities.append(V[row][column])
        if count_max(probabilities):
            katastaseis.append(Q[multiple_max_possibility(probabilities,trans,Q.index(katastaseis[column-1]),column)])
        else:
            katastaseis.append(Q[probabilities.index(max(probabilities))])
    return katastaseis,10**(max(probabilities))


priors = np.array([0.5, 0.5])
transition = np.array([[0.9, 0.1],
                        [0.1, 0.9]])
emission = np.array([[0.4, 0.4, 0.1,0.1],
                     [0.2, 0.2, 0.3,0.3]])

print("The HMM with array representation:")
print_HMM(priors,transition,emission)

#Μετατροπή πιθανοτήτων σε λογαριθμικές με βάση το 10
priors = np.log10(priors)
transition = np.log10(transition)
emission = np.log10(emission)

print("\nLogarithmic probabilities with base 10\n")
print_HMM(priors,transition,emission)

seq = "GGCT"
s,p = viderbi_algorithm(priors,transition,emission,sequence_str_to_numeric_list(seq))
print("The most probable path is:",states_list_to_str(s))
print("Its probability is (we used log10):",p)