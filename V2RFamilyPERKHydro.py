from Bio import Entrez, SeqIO
import pylab
import pandas as pd
from pylab import *

Entrez.email="moralwintertiger@gmail.com"
#attach an IdList to handle:
handle = Entrez.esearch(db="protein", term="Vmn2*[gene] NOT partial AND mus[ORGN]", retmax='250')
record = Entrez.read(handle)
#handle.close()
id_list =record['IdList']
num_residues = 900

#use handle2 to fetch sequences:
handle2 = Entrez.efetch(db="protein", id=id_list,rettype="gb", retmode="text")

#make a list of the sequences in the handle2:
input_list = []
for seq_record in SeqIO.parse(handle2,"gb"):
    input_list.append(str(seq_record.seq))

kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
       'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
       'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
       'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

kP = { 'A': 7.1,'R':5,'N':2.7,'D':2,'C': 2,
       'Q':2,'E':2.2,'G':4.4,'H':2.2,'I': 9.1,
       'L': 11.3,'K':5.2,'M': 4.1,'F': 5.7,'P':3.9,
       'S':9.4,'T':7.2,'W':1.1,'Y':3.5,'V': 9.7 }

#the key here is to have an embedded list that is being rewritten:

output_list = []
for seq in input_list:
    placeholder = []
    for letter in seq:
        placeholder.append(kd[letter])
    output_list.append(placeholder)

output_listA = []
for seq in input_list:
    placeholder = []
    for letter in seq:
        placeholder.append(kP[letter])
    output_listA.append(placeholder)

#create dictionary
'''each position in the list is a key, and the value is all of its values across lists'''
position = {}
for li in output_list:
    for index in range(0,len(li)):
        if index in position:
            position[index].append(li[index])
        else:
            position[index] = [li[index]]

positionA = {}
for li in output_listA:
    for index in range(0,len(li)):
        if index in positionA:
            positionA[index].append(li[index])
        else:
            positionA[index] = [li[index]]

#take average at each index point and add to new list of means
mean_list = []
for key in position:
    avg = mean(position[key])
    mean_list.append(avg)

mean_listA = []
for key in positionA:
    avg = mean(positionA[key])
    mean_listA.append(avg)

new_val = pd.DataFrame(mean_list)
new_valA = ((pd.DataFrame(mean_listA)))

window6 = new_val.rolling(21, center=True, win_type="triang").mean()

window5 = new_valA.rolling(21, center=True).mean()
window5 = window5*.33

x_data = range(1, len(window6)+1)



pylab.plot(x_data, window6, linewidth=2.0, label="Hydrophobicity", color="green")
pylab.plot(x_data, window5, linewidth=2.0, label="P-binding", color="blue")
pylab.legend(loc='upper left')
axis(xmin = 1, xmax = num_residues)
xlabel("Residue Number")
ylabel("Hydrophobicity")
title("Hydrophobicity and predicted P-binding for V2Rs")


handle.close()
print(len(id_list))
