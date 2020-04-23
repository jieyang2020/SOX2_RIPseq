
import sys
import string

fin = sys.stdin
fout = sys.stdout
#ferr = sys.stderr

def RevComp(seq):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = seq.translate(complements)[::-1]
	return rcseq

# iterate alignment records/lines.
umi0 = ""
reads = []
for line in fin:
	if line[0] == "@" : # if header part, skip.
		fout.write(line)
		continue

	line = line.strip().split("\t")

	# extract umi and umiQual form qname.
	temp = line[0].split(":")
	umiQual = temp[0]
	umi = temp[1]

	# restore the original qname.
	line[0] = ":".join(temp[3:])[1:] # orginal qname without '@'

	# skip if UMI has lower sequence quality.
	if umiQual == "FAIL" :
		continue

	# start main processing
	if umi0 == "" : # initialization
		umi0 = umi
		reads.append(line)
		continue
	elif umi == umi0 : # collecting the reads with same UMIs
		reads.append(line)
		continue
	else: # before starting with new UMI, process the reads with the same UMIs.
		deduppedIdx = {}
		i = 0
		while i < (len(reads)-1) :
			# if qnames are different, go to the next record.
			# Note that this program discard orphan read that is uniquely mapped.
			if reads[i][0] != reads[i+1][0]:
				i = i + 1
				continue

			# if chromosomes are different, go to the next pair of record.
			# Note that this program discard chimera.
			if reads[i][2] != reads[i+1][2]:
				i = i + 2
				continue

			# define position for the read pair
			chrom = reads[i][2]
			coord = min(int(reads[i][3]), int(reads[i+1][3]))
			pos = chrom + ":" + str(coord)

			# sequence for the read pair
			flag_i = int(reads[i][1])
			flag_ii = int(reads[i+1][1])
			if flag_i & 0b1000000 == 0b1000000: # i is first mate, i+1 is second mate
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
			if pos in deduppedIdx.keys(): # if there is a read pair with the same position,
				if len([d for d in deduppedIdx[pos] if d[1] != seq]) == len(deduppedIdx[pos]): # add only when all the existing read pairs are different from the current read pair.
					deduppedIdx[pos].append([i, seq])
			else:
				deduppedIdx[pos] = [[i, seq]]

			# move to the next pair.
			i = i + 2

		# print out the dedupped read pair.
		for idxList in deduppedIdx.values():
			for idx in idxList:
				fout.write("\t".join(reads[idx[0]])+"\n") # mate 1
				fout.write("\t".join(reads[idx[0]+1])+"\n") # mate 2

		# after the current UMI processing is done, update umi.
		umi0 = umi
		reads = [line]

fin.close()
fout.close()
