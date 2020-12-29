import sys

docstring='''Given two files (a features.bed file and a reads.bed file), append the strings from the last column of the features file to an array within an array containing information in the reads file. These strings should append only if the intervals in the second and third columns of the features file overlap the intervals of the first and second columns of the reads file.
USAGE
bederlap_finder.py <features.bed> <reads.bed>

--------------EXAMPLE-----------------
reads.bed:
chromA  10      69      read1
chromA  10      35      read2
chromA  10      55      read3
chromA  15      69      read4
chromA  80      119     read5
chromA  80      111     read6
chromA  90      119     read7
chromA  101     119     read8

feats.bed:
chromA  10      19      feat1
chromA  30      39      feat2
chromA  50      69      feat3
chromA  80      89      feat4
chromA  100     119     feat5

CALL
bederlap_finder.py features.bed reads.bed

OUTPUT
[10, 69, 'read1', ['feat1', 'feat2', 'feat3']]
[10, 35, 'read2', ['feat1', 'feat2']]
[10, 55, 'read3', ['feat1', 'feat2', 'feat3']]
[15, 69, 'read4', ['feat1', 'feat2', 'feat3']]
[80, 119, 'read5', ['feat4', 'feat5']]
[80, 111, 'read6', ['feat4', 'feat5']]
[90, 119, 'read7', ['feat5']]
[101, 119, 'read8', ['feat5']]
'''

if len(sys.argv) !=3:
        sys.exit(docstring)

feat_bed=open(sys.argv[1])
read_bed=open(sys.argv[2])

#feat_bed=open("feats.bed","r")
#read_bed=open("reads.bed","r")


feat_lines=feat_bed.readlines()

for line in read_bed.readlines():
        line=line.strip()
        line=line.split("\t")
        read=[int(line[1]),int(line[2]),str(line[3]),[],[],[]]
        for feat in feat_lines:
                feat=feat.strip()
                feat=feat.split("\t")
                if int(read[0]) <= int(feat[1]) <= int(read[1]) and int(feat[1]) <= int(feat[2]) <= int(read[1]):
                        read[3].append(str(feat[3]))
                        bpoverlap = int(feat[2]) - int(feat[1]) + 1
                        read[4].append(bpoverlap)
                        frac = float(bpoverlap) / (int(feat[2]) - int(feat[1]) + 1)
                        read[5].append(frac)
                elif int(read[0]) <= int(feat[1]) <= int(read[1]) and int(feat[2]) >= int(read[1]):
                        read[3].append(str(feat[3]))
                        bpoverlap = int(read[1]) - int(feat[1]) + 1
                        read[4].append(bpoverlap)
                        frac = float(bpoverlap) / (int(feat[2]) - int(feat[1]) + 1)
                        read[5].append(frac)
                elif int(feat[1]) <= int(read[0]) and int(read[0]) <= int(feat[2]) <= int(read[1]):
                        read[3].append(str(feat[3]))
                        bpoverlap = int(feat[2]) - int(read[0]) + 1
                        read[4].append(bpoverlap)
                        frac = float(bpoverlap) / (int(feat[2]) - int(feat[1]) + 1)
                        read[5].append(frac)
                elif int(feat[1]) <= int(read[0]) and int(feat[2]) >= int(read[1]):
                        read[3].append(str(feat[3]))
                        bpoverlap = int(read[1]) - int(read[0]) + 1
                        read[4].append(bpoverlap)
                        frac = float(bpoverlap) / (int(feat[2]) - int(feat[1]) + 1)
                        read[5].append(frac)
        print read, 'features =',len(read[3])
