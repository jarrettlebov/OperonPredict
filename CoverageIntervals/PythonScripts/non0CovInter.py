import sys

#open file
depth_file=open(sys.argv[1])

#establish an interable variable that stores each line of the opened file
pos=depth_file.readlines()

#initiate start and stop arrays to hold start and stop coordinates for each coverage interval
starts=[]
ends=[]

#add the first position in the file to the 'starts' array
starts.append(pos[0].strip("\n"))

#look for examples in the file where a given line is > than the value of the previous line+1, and append lines to 'starts' and previous lines to 'ends'
prev=next(iter(pos))
for line in pos:
	line=line.strip("\n")
	if int(line) > int(prev)+1:
		starts.append(line)
		ends.append(prev)
	prev = line

#add last position in the file to the 'ends' array
ends.append(pos[len(pos)-1].strip("\n"))

#print each array side by side with a tab separating each array
print("starts\tends")
for i in range(len(starts)):
	print(starts[i] + '\t' + ends[i])
