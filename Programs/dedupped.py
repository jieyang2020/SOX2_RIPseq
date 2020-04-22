
import sys
import string

fin = sys.stdin
fout = sys.stdout
#ferr = sys.stderr

def RevComp(seq):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = seq.translate(complements)[::-1]
	return rcseq

umi0 = ""
reads = []
for line in fin:
	if line[0] == "@" :
		fout.write(line)
		continue

	line = line.strip().split("\t")

	# umi, umiQual, original line
	temp = line[0].split(":")
	umiQual = temp[0]
	umi = temp[1]
	line[0] = ":".join(temp[3:])[1:] # don't forget to remove @

	if umiQual == "FAIL" :
		continue

	if umi0 == "" : # initialization
		umi0 = umi
		reads.append(line)
		continue
	elif umi == umi0 : # collecting the reads with same UMIs
		reads.append(line)
		continue
	else: # after collecting the reads with same UMIs,
		deduppedIdx = {}
		for i in range(0, len(reads), 2): # iterate two consecutive reads (mates): i and i+1
			
			# check qname
			if reads[i][0] != reads[i+1][0]:
				continue

			# check chromosome 
			if reads[i][2] != reads[i+1][2]:
				continue
			
			# position
			chrom = reads[i][2]
			coord = min(int(reads[i][3]), int(reads[i+1][3]))
			pos = chrom + ":" + str(coord)

			# sequence: read1 + read2
			flag_i = int(reads[i][1])
			flag_ii = int(reads[i+1][1])
			if flag_i & 0b1000000 == 0b1000000: # i: first mate, i+1: second mate
				seq = reads[i][9]
				if flag_i & 0b10000 == 0b10000: # reverse strand
					seq = RevComp(seq)
				if flag_ii & 0b10000 == 0b10000: # is second mate reverse strand just in case...
					seq = seq + RevComp(reads[i+1][9])
				else:
					seq = seq + reads[i+1][9]
			else: # i: second mate, i+1: first mate
				seq = reads[i+1][9]
				if flag_ii & 0b10000 == 0b10000: # reverse strand
					seq = RevComp(seq)
				if flag_i & 0b10000 == 0b10000: # is first mate reverse strand just in case...
					seq = seq + RevComp(reads[i][9])
				else:
					seq = seq + reads[i][9]
			
			# checking duplicates
			if pos in deduppedIdx.keys(): # update
				if len([d for d in deduppedIdx[pos] if d[1] != seq]) == len(deduppedIdx[pos]):
					deduppedIdx[pos].append([i, seq])
				#else:
				#	ferr.write(reads[i][0]+"\n")
			else: # initialization and adding the reads with same position
				deduppedIdx[pos] = [[i, seq]]

		# print out
		for idxList in deduppedIdx.values():
			for idx in idxList:
				fout.write("\t".join(reads[idx[0]])+"\n") # mate 1
				fout.write("\t".join(reads[idx[0]+1])+"\n") # mate 2

		# update
		umi0 = umi
		reads = [line]

fin.close()
fout.close()









