#!/usr/bin/env python

inter_bed=open("f_K12_intergene.bed","r")
nomRNA_bed=open("K12_non-t-rRNA.bed","r")

nomRNA=[] #set up empty nomRNA array
for r in nomRNA_bed.readlines(): # open for loop which will iterate through the lines in the nomRNA_bed file
        r=r.strip() #strip new line characters
        r=r.split("\t") #split elements in each line by \t
        nomRNA.append(int(r[1])) #append elements in the second column (first column = [0]) to empty array nomRNA[]

excludes=[]
inters=[]
for i in inter_bed.readlines():
        i=i.strip()
        i=i.split("\t")
        inters.append("#".join(i))
        for element in nomRNA:
                if int(i[1]) < element < int(i[2]):
                        #print("\t".join(i))
                        #print element
                        if "#".join(i) not in excludes:
                                excludes.append("#".join(i))

real = set(inters) - set(excludes)
print("\n".join(real))
