#!/usr/local/bin/python

##
#
# 08/21,02
#  Implement compare_lists()
# 09/23,02
#  Problem with DOS style files. In get_sequences, "\r\n" is explicitly detected
#  but should use a generate function to convert file first.
# 01/02,03
#  The function get_group_seq don't output sequences in the order of the input
#  list. Make it so.
# 05/27,03
#  This is very significant so I am writing this down. A new code seqment for
#  fasta_to_dict is written, now it proceed at BLAZING SPEED.
##

import os, sys, FileUtility, string, Translation

class fasta_manager:

	def __init__(self):
		pass

	#
	# @param coords  [seq_id][L][R]
	#
	def mask(self,fasta,coords):
		
		print ("Read fasta into dict...")
		F = self.fasta_to_dict(fasta,dflag=0)
		
		print ("Read coords into dict...")
		C   = {}
		inp = open(coords)
		inl = inp.readline()
		while inl != "":
			L = inl.split("\t")
			if C.has_key(L[0]):
				C[L[0]].append([int(L[1]),int(L[2])])
			else:
				C[L[0]] = [[int(L[1]),int(L[2])]]
			inl = inp.readline()
		
		print ("Concatenate seq...")
		for i in C:
			if F.has_key(i):
				print ("",i,"%i features" % len(C[i]))
				c = 0
				for j in C[i]:
					if c % 100 == 0:
						print (" %i x100" % (c/100))
					c += 1
					F[i] = "%s%s%s" % \
					(F[i][:j[0]-1],"N"*(j[1]-j[0]+1),F[i][j[1]:])
			else:
				print ("SEQ ABSENT:",i)
		
		print ("Write sequences...")
		self.dict_to_fasta(F,fasta+".mask")
		print ("Done!")
		
		
	def dict_to_fasta(self,fdict,fname):
		
		oup = open(fname,"w")
		for i in fdict:
			oup.write(">%s\n" % i)
			s = fdict[i]
			c = 0
			while c < len(s):
				oup.write("%s\n" % s[c:c+80])
				c += 80
		oup.close()


	##
	# This is written to get stretch from a LARGE sequence file with a single
	# sequence, such as the chromosome sequence. The other method takes forever
	# because the fasta is concatenated into a single string. A SERIOUS negative
	# for this function is that the coords CANNOT overlap.
	#
	# @param fasta  the implementation only work for one id one seq in the file
	# @param coords separated by ",". Sequence obtained is inclusive of the
	#               coords OR a file name can be passed. The file contains
	#               a pair of coords separated by "," in each line
	# @param call   whether return the sequence as dict [1] or written [0]
	# @param isfile coords in a file
	##
	def get_stretch2(self,fasta,coords,call=0,isfile=0):
		
		print ("Fasta   :",fasta)
		print ("Coord   :",coords)
		print ("isfile  :",isfile)

		# check if coord is a file
		if isfile:
			inp = open(coords,"r")
			inl = inp.readline()
			
			# first read to a dict to sort the coords
			order = {}
			while inl != "":
				L = self.rmlb(inl).split(",")
				# assume coord are not overlapping
				order[int(L[0])] = int(L[1])
				inl = inp.readline()
			
			okeys = order.keys()
			okeys.sort()
			clist = []
			for i in okeys:
				clist.extend([i,order[i]])
			inp.close()
		else:
			clist = coords.split(",")
			clist = [int(clist[0]),int(clist[1])]
		
		inp = open(fasta,"r")
		idx = self.rmlb(inp.readline())[1:]
		
		# rid of anything after the 1st space
		if idx.find(" ") != -1:
			idx = idx[:idx.find(" ")]
		
		inl = inp.readline()
		c     = 0
		get   = 0
		odict = {}
		ostr  = ""
		end   = 0
		read  = 1
		
		saved = 0
		while inl != "":
			inl = self.rmlb(inl)
			#print [inl]
			
			# only increment if this line has been read
			if read:
				c += len(inl)
			else:
				read = 1
			
			if int(clist[0]) < int(clist[1]):
				cL = int(clist[0])
				cR = int(clist[1])
			else:
				cR = int(clist[0])
				cL = int(clist[1])
			
			if c > cL:				
				if get == 0:
					#print "%s-%s" % (clist[0],clist[1]),cL,cR
					#print " start:",
					if c < cR:
						ostr += "%s\n" % inl[cL-c+len(inl)-1:]
						get = 1
					else:
						#print "end1"
						ostr += "%s\n" % inl[cL-c+len(inl)-1:\
											 cR-c+len(inl)]
						end = 1
				else:
					if c < cR:
						#print ">%s\n" % inl
						#print ".",
						ostr += "%s\n" % inl
					else:
						#print "end2"
						ostr += "%s\n" % inl[:cR-c+len(inl)]
						end = 1
										
				# store seq into dict and reset everything
				if end:
					#print " reset"
					#print " store->",[ostr]
					
					if saved % 10000 == 0:
						print (" %i x 10k" % (saved/10000))
					saved += 1
					
					# deal with coords in reverse
					if int(clist[0]) > int(clist[1]):
						ostr = trans.rc(ostr)
						print ("rc:",[ostr])
					
					odict["%s|%s" % (clist[0],clist[1])] = ostr					
					# reset everything
					end    = 0
					get    = 0
					ostr   = ""
					clist = clist[2:]
					if clist == []:
						break
					# just in case that the next entry is also in this line.
					if c > cL:
						#print ">>>no_read"
						read = 0
			# read line or not						
			if read:		
				#print " new_line"
				inl = inp.readline()		
		
		if call:
			return odict
		else:
			if isfile:
				oup = open(coords+".seg.fa","w")
			else:
				oup = open(fasta +".seg.fa","w")
			
			for i in odict:
				oup.write(">%s_%s\n%s" % (idx,i,odict[i]))
				
			print ("Done!")	
	
	# for multiple chr, derived from get_stretch3.
	# 2/24/12, coordinates used to be [chr][L][R] or can be multiple coordinates
	#          separated by ",". now allow only [seq][L1,R1,...]
	#
	# @param fasta  can have multiple sequences in one file
	# @param coords [seq][L1,R1,L2,R2...]. If L > R, take reverse complement.
	# @param seqid  Use 4th column of the coords file as sequence ID if available.
	#               (1) or not (0, default)
	# use get_stretch4 where you have Chr\tseq1,seq2\tname
	def get_stretch4(self,fasta,coords,seqid):
		
		print ("Sequence to dict...")
		seq = self.fasta_to_dict(fasta,0)
		# fasta_to_dict got rid of "\n" already
	
		print (seq.keys())
		new_dict={}
		for i in seq.keys():
		    newseq= i.strip().split(" ")[0]
		    print (newseq)
		    value= seq[i]
		    new_dict[newseq]=value
		
	
		c = 0 # count total
		m = 0 # count not in fasta
		if coords.find(",") == -1:
			print ("Read coordinates...")
			inp = open(coords)
			oup = open(coords+".fa","w")
			oup2= open(coords+".missing","w")
			inl = inp.readline()
			while inl != "":    # Go through each coord
				if c % 1000 == 0:
					print (" %i k" % (c/1000))
				c  += 1
				L  = inl.strip().split("\t")
				seqName = L[0]	# Sequnece name
				if seqName in new_dict:
					if len(L) >= 2:
						# L = [name, L, R], some may have 4th col which is IDs to be given.
						
						if L[1].find(",") == -1:
							# Deal with reverse ori
							ori= 1; cL = int(L[1]); cR = int(L[2])
							if cL > cR:
								ori = -1
							# Get sequence
							if ori == -1:
								S = new_dict[seqName][cR-1:cL]
							else:
								S = new_dict[seqName][cL-1:cR]
	
							if S == "":
								print ("ERR COORD: %s,[%i,%i]" % (seqName,cL,cR))
							else:
								if ori == -1:
									S = trans.rc(S)
								
								# If there is 4th column, use them as sequence IDs.
								if len(L) == 4 and seqid:
									oup.write(">%s\n%s\n" % (L[3],S))
								else:
									oup.write(">%s|%i-%i\n%s\n" %	(seqName,cL,cR,S))
						# name <\t> "L1,R1,L2,R2..." <\t> whatever
						elif L[1].find(",") != -1:
							genename= L[2]
							coordList = L[1].split(",")
							S = ""
							# Set orientation, only consider the first pair
							if int(coordList[0]) < int(coordList[1]):
								ori = 1
							else:
								ori = -1
							for j in range(0,len(coordList),2):
								cL = int(coordList[j]); cR = int(coordList[j+1])
								if ori == -1:
									cL = int(coordList[j+1])
									cR = int(coordList[j])
									S = new_dict[seqName][cL-1:cR] + S
								else:
									S += new_dict[seqName][cL-1:cR]
							oup.write(">%s %s|%s\n%s\n" % (genename, seqName,"-".join(coordList),S))
						else:
							print ("Unknown cooord format:",L)
							print ("Quit!")
							sys.exit(0)

				else:
					m += 1
					oup2.write(inl)
				inl = inp.readline()
			print ("Total coords:",c)
			print ("Not in seq  :",m)
			oup2.close()

		# coordinates are passed
		else:
			coords = coords.split(",")
			print ("Coords:",coords)
			oup = open("%s_%s.fa" % (fasta,"-".join(coords)),"w")
			C = []
			for i in coords:
				C.append(int(i))
			for i in seq:
				s = ""
				for j in range(0,len(C),2):
					s += seq[i][C[j]-1:C[j+1]]
				oup.write(">%s\n%s\n" % (i,s))
					
		oup.close()
		print ("Done!")	
		
	#
	# Well... get_strech2 is problematic when overlapping sequences are needed.
	# Now this method is much smarter... and much shorter too.
	# 07/19,05. New lines are not formed in some entries. Fixed?
	# 04/27,06. Reverse complement coordinates are not reflected. Fixed.
	#           Modify this for multiple sequences
	#
	# @para fasta  file with single sequence
	# @para coords coord file, tab-delim, the last two token have to be coords.
	#
	def get_stretch3(self,fasta,coords):
		
		print ("Covert fasta into a string...")
		inp = open(fasta)
		inl = inp.readlines()
		idx = self.rmlb(inl[0])[1:]
		seq = string.join(inl[1:],"")

		# get rid of line breaks, sometimes \n and \r\n are mixed.
		if seq.find("\r\n") != -1:
			seq = string.join(seq.split("\r\n"),"")
		if seq.find("\n") != -1:
			seq = string.join(seq.split("\n"),"")
		print (len(seq))
		
		print ("Read coords and output seq...")
		inp = open(coords)
		oup = open(coords+".fa","w")
		inl = inp.readline()
		c = 0
		while inl != "":
			print ([inl])
			if c % 10000 == 0:
				print (" %i x 10k" % (c/10000))
			c += 1
			L = inl.split("\t")
			
			# deal with reverse ori
			cL = int(L[-2])
			cR = int(L[-1])
			ori = 1
			if int(L[-2]) > int(L[-1]):
				cL = int(L[-1])
				cR = int(L[-2])
				ori = -1
			
			S = seq[cL-1:cR]
			if S == "":
				print ("ERR COORD: [%i,%i]" % (cL,cR))
			else:
				if ori == 1:
					oup.write(">%s|%i-%i\n%s\n" % (idx,cL,cR,S))
				else:
					S = trans.rc(S)
					oup.write(">%s|%i-%i\n%s\n" % (idx,cR,cL,S))
					#print ">%s|%i-%i %s" % (idx,cR,cL,S[:10])
						
			inl = inp.readline()
		
		print ("Done!")
		

	##
	# coord format: from id<space>L-R to id_((L-1)*3+1)_R
	## 
	def convert_header(self,fasta):

		inp = open(fasta,"r")
		oup = open(fasta+".coord.fa","w")
		inl = inp.readline()
		while inl != "":
			if inl[0] == ">":
				if inl.find(" ") == -1 or inl.find("-") == -1:
					print ("Wrong descriptor format, Quit!")
					sys.exit(0)
					
				L = inl.split(" ")
				left = int(L[1][:L[1].find("-")])
				righ = int(L[1][L[1].find("-")+1:-1])

				left = (left-1)*3+1
				righ = righ*3
				oup.write("%s_%i_%i\n" % (L[0],left,righ))

			else:
				oup.write(inl)
			inl = inp.readline()

		
	# rid of coords
	#
	# @style  coord connected to id via underscore: 0 [default], or space [1]
	def delete_coord(self,fasta,style):

		print ("NOT SURE THIS IS WORKING YET")

		inp = open(fasta,"r")
		oup = open(fasta+".mod.fa","w")

		inl = inp.readline()
		while inl != "":
			if inl[0] == ">":
				if style == 0:
					if inl.find("_") == -1:
						oup.write(inl)
				elif style == 1:
					if inl.find(" ") == -1:
						oup.write(inl)


			inl = inp.readline()


	##
	# Return a dict with seq id as key, seq as value. Will convert things to
	# UNIX format without "\r".
	#
	# 11/15 Major rewrite of the method
	# 11/17 Deal with redundant id, before they are simply replaced.
	#
	# @param fasta
	# @param desc    include desc[1] or not[0,default]. The description is
	#	             separated from seq id with a space. If include, it will
	#	             be the first element of the list, the second will be seq
	# @param verbose display line count [1], or not [0, default]
	# @param dflag   delimit by space or not, default no
	# @param newline rid of new line or not
	##
	def fasta_to_dict(self,fasta,dflag=0,verbose=0,newline=0):
	
		# each idx should only occur once.
		inp    = open(fasta,"r")
		inl = inp.readline()
		fdict  = {} # idx as key, seq as value
		c = 0
		N = 0 
		while inl != "":
			inl = self.rmlb(inl)
			next = 0
			idx  = ""
			desc = ""
			if inl == "":
				pass
			elif inl[0] == ">":
				if verbose and c%1e3 == 0:
					print (" %i k" % (c/1e3))
				c += 1
				# rid of anything after space if asked
				if dflag and inl.find(" ") != -1:
					desc = idx[idx.find(" ")+1:]
					idx = idx[:idx.find(" ")]
				else:
					idx = inl[1:]
				# count lines and store seq into a list
				slist = []
				inl = inp.readline()
				while inl[0] != ">":
					inl = inl.strip()
					# Add new line char, do this after strip because I do not
					# want to have cases of /r/n 
					if newline:
						inl = inl + "\n"
					slist.append(inl)
					inl = inp.readline()
					if inl == "":
						break
				seq = "".join(slist)
				
				if idx in fdict.keys():
				#if fdict.has_key(idx):
					if verbose:
						print ("Redundant_id:",idx,)
					if dflag:
						if len(fdict[idx][1]) < len(seq):
							fdict[idx] = [desc,seq]
							if verbose:
								print ("longer")
						else:
							if verbose:
								print ("shorter")
					else:
						if len(fdict[idx]) < len(seq):
							fdict[idx] = seq
							if verbose:
								print ("longer")
						else:
							if verbose:
								print ("shorter")
				else:
					N += 1
					if dflag:
						fdict[idx] = [desc,seq]
					else:
						fdict[idx] = seq
				next = 1
			
			# so no extra line is read, because of the innder while
			if not next:
				inl = inp.readline()
		inp.close()
		if verbose:
			print ("Total %i sequences, %i with non-redun names" % (c,N))
		
		#print(fdict)
		return fdict

	#
	# @param t  Truncated or not, default no.
	def fasta_to_phylip(self,fasta,t=0):
		
		print ("To oneline...")
		self.fasta_to_oneline(fasta)
		
		print ("To phylip...")
		inp = open(fasta+".pep")
		oup = open(fasta+".phylip","w")
		inl = inp.readlines()
		c   = 0
		
		# Figure out the longest name
		nlen = 0
		for i in inl:
			n = i.split("\t")[0]
			if len(n) > nlen:
				nlen = len(n)
				#print nlen,n
		#print nlen

		for i in inl:
			L = i.strip().split("\t")
			
			# Write header
			if c == 0:
				oup.write("%i %i\n" % (len(inl),len(L[1])))
				c += 1
			n = L[0]
			if t:
				if len(n) > 10:
					print ("Truncate name:",n)
					n = n[:10]
				oup.write("%s%s%s\n" % (n," "*(10-len(n)+1),L[1]))
			else:
				oup.write("%s%s%s\n" % (n," "*(nlen-len(n)+1),L[1]))
			
		inp.close()
		oup.close()
		print ("Done!")

	#
	# This is mainly done for HMMER3...
	#
	def fasta_to_stockholm(self,fasta):
		fa = self.fasta_to_dict(fasta)
		blocklen = 80
		
		# Find out length of the longest name and number of blocks
		nlen = 0
		nblc = 0
		for i in fa:
			if len(i) > nlen:
				nlen = len(i)
			if nblc == 0:
				nblc = len(fa[i])/80.0
				if nblc > int(nblc):
					nblc = int(nblc)+1
				else:
					nblc = int(nblc)
				
		# Output sequence blocks	
		oup = open(fasta+".stockholm","w")
		oup.write("# STOCKHOLM 1.0\n\n")
		
		nkeys = fa.keys()
		nkeys.sort()
		for i in range(nblc):
			for j in nkeys:
				n = "_".join(j.split(" ")) # Rid of space in names
				#print nlen, len(n)
				oup.write("%s%s%s\n" % (n,
				                        " "*(nlen-len(n)+2),
				                        fa[j][blocklen*i:blocklen*(i+1)]))
			oup.write("\n")
			
		oup.write("//")
		oup.close()


	##
	# Generate fasta file based on species passed, assuming desc like:
	#  >XX_id
	#
	# @param fasta
	##
	def get_sp(self,fasta,species):

		print ("Get sequences for a particular species:")
		print (" Fasta  :",fasta)
		print (" Species:",species)
		
		inp = open(fasta,"r")
		oup = open(fasta+".%s.fa" % species,"w")
		inl = inp.readline()
		slist = []

		w = 1 # write or not
		c = 0
		while inl != "":
			if inl[0] == ">":
				if inl[1:3] == species:
					w = 1
					c += 1
				else:
					w = 0

			if w:
				oup.write(inl)       
			inl = inp.readline()

		print ("Total %i sequences" % c)

	##
	# Count the number of sequences for species specified
	#
	# @fasta
	# @sp     Species header file, each line one species.
	##
	def count_sp(self,fasta,sp):
		print ("Read species info into a dict...")
		inp = open(sp)
		inl = inp.readlines()
		S   = {}
		for i in inl:
			S[i.strip()] = 0
		
		print ("Go through fasta...")
		inp = open(fasta)
		inl = inp.readline()
		while inl != "":
			if inl[0] == ">":
				n  = inl.strip()[1:]
				sp = []
				for j in S:
					if j in n:
						sp.append(j)
				if len(sp) == 0:
					print (" No sp definition:",n)
				elif len(sp) >1:
					print (" Ambiguous",n,sp)
				else:
					S[sp[0]] += 1
			
			inl = inp.readline()
		
		print ("Generate output...")
		oup = open(fasta+".sp_count","w")
		for i in S:
			oup.write("%s\t%s\n" % (i,S[i]))
		
		print ("Done!")
	
	#
	# Convert GFF to coord file
	#
	def gff_to_coord(self,gff):
		inp = open(gff)
		oup = open(gff+".coord","w")
		inl = inp.readline()
		while inl != "":
			#print (inl)
			T = inl.strip().split("\t") # tokens
			if len(T) > 6:
			    #print (T)
			    C = T[0]					# chr
			    if T[6] == "+":				# orientation
			        L = T[3]				# left coord
			        R = T[4]				# right coord
			    else:
			        L = T[4]
			        R = T[3]
			    N = ""						# sequence name
			    n = T[-1].split(";")
			    if "Name" in T[-1]:			# has Name tag, # refers to last item in list
			        for j in n:
			            if "Name" in j:
			                N = j.split("=")[1]
			                break
			            elif T[-1] != "":			# no name tag but not empty, use 1st
			                N = n[0].split("=")[1]
			            else:
			                print ("No desc:",T)
			    if N != "":
			        oup.write("%s\t%s,%s\t%s\n" % (C,L,R,N))
			
			else:
			    print(T, "line not used")
			inl = inp.readline()
			
		print ("Done!")
#
	# Convert GFF to coord file for 1000bp upstream of genes
	#
	def gff_prom_to_coord2(self,gff):
		inp = open(gff)
		oup = open(gff+"_prom.coord","w")
		inl = inp.readline()
		tmp_list=[]
		count= 0
		while inl != "":
			#print (inl)
			T = inl.strip().split("\t") # tokens
			if len(T) > 6:
			    C = T[0]					# chr
			    #print (T[2])
			    if T[2] == 'mRNA' and count== 0:
			        #count= count+1
			        tmp_list=[]
			        if T[6] == "+":				# orientation
			            L = T[3]				# left coord
			            R = T[4]				# right coord
			            prom= int(L)-1000                   #get 1000 bp upstream of L
			        else:
			            L = T[4]
			            R = T[3]
			            prom = int(L)+1000
			        #print (prom, L)
			        if prom > 0:
			            tmp_list=[L, prom]
			        else:
			            tmp_list=[L, 0]
			    #if T[2]=='mRNA' and count== 1:            #only get the protein name
			        N = ""				     # sequence name
			        n = T[-1].split(";")     # get gene name
			        #print(n)
			        if "Name" in T[-1]:			# has Name tag, # refers to last item in list
			            for j in n:
			                if "Name" in j:
			                    N = j.split("=")[1]
			                    #print(N)
			                    break
			                elif T[-1] != "":			# no name tag but not empty, use 1st
			                    N = n[0].split("=")[1]
			                else:
			                    print ("No desc:",T)
			        if N != "":
			            #print(tmp_list)
			            L2= tmp_list[0]
			            prom2= tmp_list[1]
			            oup.write("%s\t%s,%s\t%s\n" % (C,prom2,L2,N))
			        count = 0
			
			else:
			    pass
			    #print(T, "line not used")
			inl = inp.readline()
			
		print ("Done!")
		
	def gff_prom_to_coord_5utr(self,gff):
		inp = open(gff)
		oup = open(gff+"_prom-5utr.coord","w")
		inl = inp.readline()
		tmp_list=[]
		count= 0
		while inl != "":
			#print (inl)
			T = inl.strip().split("\t") # tokens
			if len(T) > 6:
			    C = T[0]					# chr
			    #print (T[2])
			    if T[2] == 'mRNA' and count== 0:
			        count= count+1
			        tmp_list=[]
			        if T[6] == "+":				# forward orientation
			            L = T[3]				# left coord
			            R = T[4]				# right coord
			            prom1= int(L)-1000                   #get 1000 bp upstream of TSS
			            prom2= int(L)+500                    #get 200bp downstream of TSS
			        else: # get coords if reverse orientation
			            L = T[4]
			            R = T[3]
			            prom1 = int(L)+1000
			            prom2 = int(L)-500
			        #print (prom, L)
			        if prom1 > 0:
			            tmp_list=[prom2, prom1]
			        else:
			            tmp_list=[prom2, 0]
			    #if T[2]=='mRNA' and count== 1:            #only get the protein name
			        N = ""				     # sequence name
			        n = T[-1].split(";")     # get gene name
			        #print(n)
			        if "Name" in T[-1]:			# has Name tag, # refers to last item in list
			            for j in n:
			                if "Name" in j:
			                    N = j.split("=")[1]
			                    #print(N)
			                    break
			                elif T[-1] != "":			# no name tag but not empty, use 1st
			                    N = n[0].split("=")[1]
			                else:
			                    print ("No desc:",T)
			        if N != "":
			            #print(tmp_list)
			            L2= tmp_list[0]
			            prom= tmp_list[1]
			            oup.write("%s\t%s,%s\t%s\n" % (C,prom,L2,N))
			        count = 0
			
			else:
			    pass
			    #print(T, "line not used")
			inl = inp.readline()
			
		print ("Done!")	
	#clear spaces	
	def clear_space(self,string): #returns tab-delimited
                string = string.strip()
                while "  " in string:
                    string = string.replace("  "," ")
                string = string.replace(" ","")
                return (string)
	def gff_promoter_to_coord(self,gff):
		inp = open(gff)
		oup = open(gff+".coord","w")
		inl = inp.readline()
		while inl != "":
		    #gene = ""
		    #print (gene)
		    if inl.startswith('##'):
		        T1 = inl.strip().split(" ") #split line
		        gene = T1[1]
		        #oup.write("%s\t" % (gene))
		    elif inl.startswith('#'):
		        T1 = inl.strip().split(" ") #split line
		        gene = T1[1]
		        #oup.write("%s\t" % (gene))
		    else:
		        T = inl.strip().split("\t") # tokens
		        C = T[0]					# chr
		        if T[6] == "+":				# orientation
		            L = T[3]				# left coord
		            R = T[4]				# right coord
		        else:
		            L = T[4]
		            R = T[3]
				
		        if self.clear_space(T[2]) == "promoter": 			# check if this is the promoter
		            oup.write("%s\t%s\t%s\t" % (C,L,R))
		            oup.write("%s\n" % (gene))
		        else:
		            pass
		    #print (gene)
		    	
		    inl = inp.readline()
			
		print ("Done!")
	# getr coords for cds region
	def gff_cds_to_coord(self,gff):
		D={}
		inp = open(gff)
		oup = open(gff+"_cds.coord","w")
		inl = inp.readline()
		tmp_list=[]
		count= 0
		while inl != "":
			T = inl.strip().split("\t") # tokens
			if len(T) > 6:
			    C = T[0]					# chr
			    if T[2] == 'CDS':
			        tmp_list=[]
			        if T[6] == "+":				# orientation
			            L = T[3]				# left coord
			            R = T[4]				# right coord
			        else:
			            L = T[4]
			            R = T[3]
			        tmp_list=[L, R]
			        # get alt splice
			        altspl= T[7]
			        N = ""				     # sequence name
			        n = T[-1].split(";")     # get gene name
			        #print(n)
			        if "ID" in T[-1]:			# has Name tag, # refers to last item in list
			            for j in n:
			                if "ID" in j:
			                    N = j.split("=")[1]
			                    #print(N)
			                    N= N.replace("cds-","")
			                    break
			                elif T[-1] != "":			# no name tag but not empty, use 1st
			                    N = n[0].split("=")[1]
			                else:
			                    print ("No desc:",T)
			                    N= N.remove("cds-")
			        if N != "":
			            #print(tmp_list)
			            L2= tmp_list[0]
			            R2= tmp_list[1]
			            if N not in D:
			            	D[N]= [C,L2,R2]
			            else:
			            	D[N].append(L2)
			            	D[N].append(R2)
			
			else:
			    pass
			    #print(T, "line not used")
			inl = inp.readline()
		print(D)
		for key in D:
			data= D[key]
			C=data[0]
			clist=data[1:]
			cstr=",".join(clist)
			oup.write('%s\t%s\t%s\n' % (C,cstr,key))
		#write out
		print ("Done!")
	#
	# Rename all files in dir, if they have .fa or .fasta extension.
	#
	def rename_all(self,targetDir,name):

		fList = os.listdir(targetDir)		
		for i in fList:
			if i[-3:] == ".fa" or i[-6:] == ".fasta":	
				self.rename2("%s/%s" % (targetDir,i),name,"")

		print ("rename_all done!")
		
	#
	# New rename function. But this function is not so good. It's slow and
	# only output seq if new name is found.. Should throw this away.
	#
	def rename2(self,fasta,name,ignore=""):
		
		print ("Read fasta file:",fasta)
		fdict = manager.fasta_to_dict(fasta,0)
				
		print ("Rename sequence...")
		inp = open(name,"r")
		inl = inp.readline()
		oup = open(fasta+"_rename.fa","w")
		absent = []
		ndict  = {}
		countN = 0
		countS = 0
		
		while inl != "":
			countN += 1
			inl = inl.strip().split("\t")
			if inl[1] in fdict.keys():
				countS += 1
				oup.write(">%s\n%s\n" % (inl[0],fdict[inl[1]]))

				if inl[1] in ndict.keys():
					ndict[inl[1]]+= 1
				else:
					ndict[inl[1]] = 1
			else:
				absent.append(inl[1])
			
			inl = inp.readline()

		
		oup = open(fasta + "_rename.log","w")
		oup.write("In name file but not in sequence file:\n ")
		oup.write("%s" % "\n".join(absent))
		
		oup.write("\nSequence id redundant:\n ")
		for i in ndict:
			if ndict[i] > 1:
				oup.write(" %s: %i\n" % (i,ndict[i]))
		
		print ("Fasta: %i entries, %i with new name" % (len(fdict.keys()),countS))
		print ("Name : %i entries, %i unique, %i not in fasta" % \
				(countN,len(ndict.keys()),len(absent)))
		
		print ("Done!")
		

	##
	# Rename id within fasta file
	#
	# @param fasta
	# @param name  name file in [new][old] format. Tab-delimited.
	# @param ignore Ignore whatever is after the passed char. Default nothing.
	#
	##
	def rename(self,fasta,name,ignore=""):

		# put names to a dict, old as key, new as value
		inp    = open(name,"r")
		inl = inp.readline()
		ndict  = {}
		countR = 0
		found = 0
		
		print ("Read name file...")
		oup_log.write("Redundant names:\n")
		while inl != "":
			inl = self.rmlb(inl)
			L = inl.split("\t")
			if ndict.has_key(L[1].lower()):
				print ("",L[1])
				found = 1
				countR = countR+1
			else:
				ndict[L[1].lower()] = L[0]
			inl = inp.readline()
		if not found:
			print (" none")

		# scan fasta descriptor line
		print ("Read fasta file and change names...")
		inp    = open(fasta,"r")
		oup    = open(fasta+"_rename.fa","w")
		inl = inp.readline()
		countF = 0  # found
		countN = 0  # not found
		oup_log.write("\nFasta name not found:\n")
		while inl != "":
			if inl[0] == ">":
				inl = self.rmlb(inl)
				if ignore == "" or inl.find(ignore) == -1:	
					if ndict.has_key(inl[1:].lower()):
						oup.write(">%s\n" % ndict[inl[1:].lower()])
						countF = countF+1
					# if no new name, just output old ones.
					else:
						oup.write(inl+"\n")
						oup_log.write(" %s\n" % inl[1:])
						countN = countN+1
				else:
					if ndict.has_key(inl[1:inl.find(ignore)].lower()):
						oup.write(">%s%s\n" % \
						         (ndict[inl[1:inl.find(ignore)].lower()],\
								  inl[inl.find(ignore):]))
						countF = countF+1
					# if no new name, just output old ones.
					else:
						oup.write(inl+"\n")
						oup_log.write(" %s\n" % inl[1:])
						countN = countN+1
					
			else:
				oup.write(inl)
			inl = inp.readline()
		
		oup.close()
		print ("Found %s, not found %s, %i redundant" % (countF,countN,countR))
		print ("Done!")

	##
	# Convert fasta seq id back to its original ones
	#
	# @param fasta      the fasta desc line can have just ">name\n" or something
	#		   like ">name description". The script looks for space as
	#		   delimiter in this case, so name should NOT contain any
	#		   of it.
	# @param desc_flag  should see string after space as desc [1] or not [0]
	# @param name_file  in the format [new_name][old_name]
	##
	def change_names(self,fasta,name_file,desc_flag,delim=" "):

		print ("Change ID in Fasta file")
		print (" Fasta:",fasta)
		print (" Name :",name_file)
		
		# read names into a dict
		ndict = f_util.file_to_dict(name_file,5)

		#print ndict
		
		# compare names against the fasta file
		inp = open(fasta,"r")
		oup = open(fasta+".mod.fa","w")
		inl = inp.readline()
		countF = 0
		countA = 0
		while inl != "":
			if inl[0] == ">":
				countA = countA + 1

				old  = inl[1:-1]
				desc = ""
				if desc_flag and inl.find(" ") != -1:
					old  = old[:old.find(" ")]
					desc = inl[inl.find(" ")+1:-1]
					print (old,">>",desc)
				
				if ndict.has_key(old):
					if desc != "":
						oup.write(">%s %s\n" % (ndict[old],desc))
					else:
						oup.write(">%s\n" % ndict[old])
					countF = countF+1
				else:
					oup.write(inl)
					print ("Unknwon id:",inl[1:-1])
			else:
				oup.write(inl)       

			inl = inp.readline()

		print ("Total %i sequences, %i with new names" % (countA,countF))
		print ("Done!")
		

	##
	# Convert names into indices in a fasta file
	#
	# @param start  the starting index
	##
	def index_names(self,fasta,start):

		print ("Index name for:")
		print (" Fasta:",fasta)
		print (" Start index:",start)
		inp    = open(fasta,"r")
		oup1   = open(fasta+".fa","w")
		oup2   = open(fasta+".name","w")
		inl = inp.readline()
		c      = start
		ndict  = {}
		while inl != "":
			if inl[0] == ">":
				seq_id = inl[1:-1]
				if not ndict.has_key(seq_id):
					ndict[seq_id] = c
					oup1.write(">%i\n" % c)
					oup2.write("%s\t%i\n" % (seq_id,c))
					c = c+1
				else:
					print ("Redundant:",seq_id)
			else:
				oup1.write(inl)

			inl = inp.readline()

		print (" Ending index:",c-1)
		print ("Done!")
		return ndict
	
	#
	# @param pair	A pair file for other purposes, like batch_bl2.
	# @param fasta	Sequence files that should cover all sequences specified in
	#               the pair file.
	#
	def index_pairs(self,fasta,pairs):
		print ("Read Fasta into dict...")
		fdict = self.fasta_to_dict(fasta)
		fk    = fdict.keys()
		fk.sort()
		print (" %i sequences" % len(fk))
		
		# Output indexed sequences and name/index
		oup = open(fasta+".indexed","w")
		oup2= open(fasta+".names","w")
		fkd = {}
		for i in range(len(fk)):
			n = fk[i]
			fkd[n] = i
			oup.write(">%i\n%s\n" % (i,fdict[n]))
			oup2.write("%s\t%i\n" % (n,i))
				
		print ("Read pair file...")
		inp = open(pairs)
		inl = inp.readline()
		oup = open(pairs+".indexed","w")
		while inl != "":
			L = inl.strip().split("\t")
			if L[0] not in fdict or L[1] not in fdict:
				print (" not in fasta:")
				if L[0] not in fdict:
					print (" ",L[0])
				if L[1] not in fdict:
					print (" ",L[1])
			else:
				c1 = fkd[L[0]]
				c2 = fkd[L[1]]
				oup.write("%i\t%i\n" % (c1,c2))
			inl = inp.readline()
		
		inp.close()
		oup.close()
		oup2.close()
		

	##
	# 
	# @param fasta
	# @param return_dict boolean flag for deciding if file output [0, default]
	#                    should be generated or a dict should be returned [1].
	# @param x           determine also the masked sequences. Only works for 
	#                    text output but not direct call.
	##
	def get_sizes(self,fasta,return_dict=0,x=0):

		inp = open(fasta,"r")
		inl= inp.readline()
		oup = open(fasta+".size","w")

		S  = 0
		X  = 0
		ID = ""
		sdict = {}
		T  = 0 # total characters
		while inl != "":
			if inl[-2:] == "\r\n":
				inl = inl[:-2]
			elif inl[-1] == "\n":
				inl = inl[:-1]			
			
			# just in case there is empty lines between sequences
			if inl == "":
				pass
			elif inl[0] == ">":
				# output size, reset variables
				if S != 0:
					if return_dict:
						sdict[ID] = S
					else:
						if x:
							pMasked = float(X)/float(S)
							oup.write("%s\t%i\t%i\t%f\n" %(ID,S,X,pMasked))
						else:
							oup.write("%s\t%i\n" %(ID,S))
						T += S
				ID = inl[1:]
				S  = 0
				X  = 0
			else:
				# Don't count stop
				if inl[-1] == "*":
					inl = inl[:-1]
				S = S + len(inl)
				X = X + len(inl.split("x")) - 1
				
			inl= inp.readline()
		
		# write the last entry
		T += S
		print ("Total length:",T)
		if return_dict:
			sdict[ID] = S
			return sdict
		else:
			if x:
				pMasked = float(X)/float(S)
				oup.write("%s\t%i\t%i\t%f\n" %(ID,S,X,pMasked))
			else:
				oup.write("%s\t%i\n" %(ID,S))
			print ("Done!")

	##
	# @param fasta  the fasta sequence file. Should include a sequence called
	#               OUT if the include_o is set to 1.
	# @param group  the tab-delimited file with group specification
	# @param create_d create dir for group. Default no [0]
	# @param includ_o include outgroup sequence or not. Default 0, no outgroup
	#               included.
	##
	def get_group_seq(self,fasta,group,create_d,includ_o):
		
		# tdict = {group:{gene:[rest]}}
		tdict = f_util.file_to_dict(group,7)
		fdict = self.fasta_to_dict(fasta)
		
		if includ_o and not fdict.has_key("OUT"):
			print ("No sequence called OUT but you want to include outgroup.")
			print ("QUIT!")
			sys.exit(0)
		
		print ("Get sequences for group:\n")
		keys = tdict.keys()
		keys.sort()
		# iterate group
		for i in keys:
			print (" ",i)
			fname = i			
			if(create_d):
				os.system("mkdir %s" % i)
				fname = "./%s/%s" % (i,i)
					
			oup = open(fname+".fa","w")
			for j in tdict[i]:
				if fdict.has_key(j):
					# notice that I only take the 2nd element of fdict[j].
					# the first contain some desc-like stuff, not added.
					oup.write(">%s\n%s\n" % (j,fdict[j]))
				else:	
					print ("   %s missing" % j)	
					
			if includ_o:
				oup.write(">OUT\n%s\n" % fdict["OUT"][1])
						
			oup.close()
		
		# this is for direct method call in AlnUtility.batch_tree()
		return keys
			
	#
	# A simpler get_sequences without reading seq into dict
	#
	# @param fasta
	# @param name   single column or tab delimited file
	# @oaram tokens the tokens to regarded as ID, separated by ',', default ""
	#	
	def getseq2(self,fasta,name,tokens="",verbose=1):
		
		#print "Read name file:",name
		if tokens == 0 or tokens == "":
			tokens = [0]
		else:
			tmp = tokens.split(",")
			tokens = []
			for i in tmp:
				tokens.append(int(i))
		
		inp = open(name)
		inl = inp.readline()
		ndict = {}
		countR = 0
		rlist  = []
		while inl != "":	
			L = inl.strip().split("\t")
			for j in tokens:
				if L[j] == "":
					pass
				elif L[j] not in ndict:
					ndict[L[j]] = 0
				else:
					countR += 1
					rlist.append(L[j])
			inl = inp.readline()
		countN = len(ndict.keys())
		
		#print "\nRead fasta file:"
		oup = open(name+".fa","w")
		countA = 0
		countF = 0
		for i in fasta.split(","):
			if verbose:
				print ("",i)
			inp = open(i)
			inl = inp.readline()
			flagW  = 0
			while inl != "":
				if ">" in inl:
					countA += 1
					# number found is the same as number of unique names
					if countF == countN:
						break					
					N = self.rmlb(inl)[1:]
					N = N.split(" ")[0]
					# in dict and has not been written
					if N in ndict and ndict[N] == 0:
						#print " ",N
						ndict[N] = 1
						if verbose:
							if countF % 1e2 == 0:
								print (" %i x100" % (countF/1e2))
						countF += 1
						flagW = 1
					else:
						flagW = 0
				if flagW:
					oup.write(inl)

				inl = inp.readline()
				
		oup.close()
	
		missing = []
		if countF < len(ndict.keys()):
			for i in ndict:
				if ndict[i] == 0:
					missing.append(i)
	
		if verbose:
			print (" Redun:%i, Unique:%i, Found:%i, Missed:%i" % \
				(countR,countN,countF,len(missing)))

		return missing
			
	##
	# All names are coverted to lower case when compared but the output file
	# retain the same case as the fasta file. If the fasta file sequence header
	# contain space, the sequence id is whatever after ">" and  before the space.
	#
	# @param segment  no[0,default], yes[1]. This is done based on the name fie
	#		 where coordinates are specified like:
	#		   >name [whatever] L-R
	#		 The method uses space as delimiter. The last token is
	#		 assumed to be the one with coord info.
	# @param type     fasta seqeunce is protein [0,default] or nucleotide [1]
	# @param match    exact match to source fasta [0, default], or just match
	#		 the passed name to the exact number of char after ">" in
	#		 fasta file [1]
	# @param call     0 (default): run whole module, or 1: method call from other
	#		 module.
	##
	def get_sequences(self,fasta,name_file,segment=0,type=0,match=0,call=0,
					  outname=""):
		
		print ("Fasta:",fasta)
		print ("Name :",name_file)
		
		# construct name list
		if not call:
			inp     = open(name_file,"r")
			inls = inp.readlines()
			names   = {}
			print ("Redundant names in name file:")
			flag = 0
			for i in inls:
				# check format
				if i[-2:] == "\r\n":
					i = i[:-2]
				elif i[-1] == "\n":
					i = i[:-1]
					
				# skip empty line
				if i == "":
					continue
				
				# see if anything else is encoded in the name
				if i.find(" ") != -1:
					llist = i.split(" ")
					#llist[-1] == llist[1][:-1] disabled due to format checking
				else:
					llist = [i]
					#llist = [i[:-1]], modified due to format checking

				if not names.has_key(llist[0].lower()):
					if len(llist) == 1:
						names[llist[0].lower()] = [0,""]
					else:
						names[llist[0].lower()] = [0,llist[1:]]
				else:
					print (" ",i)
					flag = 1
			
			if not flag:
				print (" None...")

			print ("Total %i names" % len(names.keys()))
		else:
			names = name_file
		
		if outname != "":
			name_file = outname
		
		countP = 0 # counting problem coord
		countO = 0 # counting ok coord
		
		# scan fasta file
		inp     = open(fasta,"r")		       
		oup     = open(name_file+".fa","w")
		inl  = inp.readline()

		count = 0
		flag = 0
		nkeys = names.keys()
		#print nkeys
		key_found = []
		while inl != "":
			found = 0
			if inl[0] == ">":
				#print inl[:-1]
				
				# convert id to lower case
				if inl.find(" ") != -1:
					id = inl[1:inl.find(" ")].lower()
				else:
					id = inl[1:-1].lower()
					if id[-1] == "\r":
						id = id[:-1]
				
				# check if any id match the same length string after ">", slow
				gotit = 0
				nname = ""
				if match:
					for i in nkeys:
						#print "",i,"vs",id[:len(i)]
						if i == id[:len(i)]:
							nname = i
							gotit = 1
							#print "  found!"
							break
				
				# include modification that won't take "redundant" sequence.
				# but this modification may impact the codes getting domains
				# need to see if this is the case.
				if (names.has_key(id) and names[id][0] != 1) or gotit:
					count = count +1
					#print "Found:",count,id
					
					if segment == 0:
						# take care of two domains in the same sequecnce
						if names[id][0] > 0:
							if inl[-2:] == "\r\n":
								inl = inl[:-2]	
							elif inl[-1] == "\n":
								inl = inl[:-1]
							oup.write("%s %i\n" % (inl,names[id][0]))						
						else:
							oup.write(inl)

						inl = inp.readline()
						while inl[0] != ">":
							oup.write(inl)
							inl = inp.readline()					 
							if inl == "":
								break

					elif segment > 0:
						coord = names[id][1][-1]
						L     = int(coord[:coord.find("-")])
						R     = int(coord[coord.find("-")+1:])
						
						if (abs(R-L)+1)%3 != 0:
							countP = countP+1
							#print "",id,L,R
						else:
							countO = countO+1

						# take care of two domains in the same sequecnce
						if names[id][0] > 0:
							if inl[-2:] == "\r\n":
								inl = inl[:-2]	
							elif inl[-1] == "\n":
								inl = inl[:-1]
							oup.write(">%s %i %s\n" % (id,names[id][0],coord))						
						else:
							oup.write(">%s %s" % (id,coord))

						inl = inp.readline()
						seq    = ""
						while inl[0] != ">":
							seq = seq + inl[:-1]
							inl = inp.readline()			 
							if inl == "":
								break					   

						if R > L:
							if type:
								seq = seq[L-1:R-3]
							else:
								seq = seq[L-1:R-1]
						else:
							if type:
								seq = seq[R-1:L-3]
							else:
								seq = seq[R-1:L-1]

						oup.write(seq+"\n")
	
					# set names value to 1
					if match:
						names[nname] = [1,names[nname][1]]
					else:
						names[id] = [1,names[id][1]]
					found = 1
					

			if not found:
				inl = inp.readline()

		# output a list of sequences not found
		oup.close()
		oup = open(name_file+".log","w")
		if not call:
			for i in names.keys():
				if names[i][0] == 0:
					flag = 1
					oup.write("%s\n" % i)

			if not flag:
				print (" None...")
			print ("Found %i out of %i" % (count,len(names.keys())))
			if type:
				print ("%i are in frame, %i out of frame" % (countO,countP))

	##
	# Convert fasta style file to one line format with [seq_id][seq]
	#
	# @param format this is only for direct method call, not from command line.
	#	       0: default, output as a file
	#	       1: output as a dict, with id as key, seq as value
	#        d delimiter, defailt the whole line will be taken.
	##
	def fasta_to_oneline(self,fasta,format=0,d=""):
		
		#print "Convert Fasta file to one line format:"
		
		inp = open(fasta,"r")
		
		inl = inp.readline()
		count = 0
		odict = {}
		seq   = ""
		id = ""
		while inl != "":
			if inl[-2:] == "\r\n":
				inl = inl[:-2]
			elif inl[-1] == "\n":
				inl = inl[:-1]
			if inl == "":
				inl = inp.readline()
				continue
		
			found = 0
			if inl[0] == ">":
				
				if seq != "":
					odict[id] = seq
					seq = ""
				
				if d != "" and inl.find(d) != -1:
					id = inl[1:inl.find(d)]
				else:
					id = inl[1:]

				odict[id] = ""
				inl = inp.readline()
				while inl[0] != ">":
					if inl[-2:] == "\r\n":
						inl = inl[:-2]
					elif inl[-1] == "\n":
						inl = inl[:-1]

					seq = seq + inl
					inl = inp.readline()
					
					if inl == "":
						break

				found = 1
				count = count +1		

			if not found:
				inl = inp.readline()

		# put the last seq into dict
		odict[id] = seq

		if format:
			return odict
		else:
			oup = open(fasta+".pep","w")
			okeys = odict.keys()
			okeys.sort()
			for i in okeys:
				oup.write("%s\t" % i)
				oup.write("%s\n"  % odict[i])
				
		print ("Total %i sequences converted" % count)


	def oneline_to_fasta(self,oneline):

		inp = open(oneline,"r")
		oup = open(oneline+".fa","w")
		inl = inp.readline()
		while inl != "":
			llist = self.rmlb(inl).split("\t")

			oup.write(">"+llist[0]+"\n")

			while len(llist[1]) > 0:
				if len(llist[1]) > 80:
					oup.write(llist[1][:80]+"\n")
				else:
					oup.write(llist[1][:80])
				try:
					llist[1] = llist[1][80:]
				except IndexError:
					break
			oup.write("\n")
			inl = inp.readline()

		print ("Output to %s.fa" % oneline)
		print ("Done!")
	
	
	##
	# Parse the description out of fasta.
	#
	# @param fasta  the fasta file with meaningful header
	# @param style  the type of fasta file:
	#		GenBank(gb) : gi|xxxxx|db|acc| desc [species]
	#		MAtDB(matdb): gene_name desc
	#		EnsEMBL(ensembl): [id] [Gene:] [Clone:] [Contig:] [Chr:] 
	#			     [Basepair:] [status:]
	##
	def parse_desc(self,fasta,style,delim):
		
		print ("Fasta:",fasta)
		print ("Style:",[style])
		print ("Delim:",[delim])
		
		inp = open(fasta,"r")
		oup = open(fasta+".desc","w")
		
		#if style == "gb":
		#       oup.write("GI#\tdatabase\tAcc\tAcc2\tDescription\tOrganism\n")
		#elif style == "ensembl":
		#       oup.write("ID\tGene\tClone\tContig\tChr\tBasePair\tStatus\n")
		
		print ("\nRead through fasta file...")
		inl = inp.readline()
		c = 0
		while inl != "":
			if inl[0] == ">":
				if c % 1000 == 0:
					print (" %i k" % (c/1000))
				c += 1
					
				# rid of ">" and "\n"
				inl = self.rmlb(inl)[1:]
				desc = species = ""

				# GenBank style
				if style == "gb":
					glist  = inl[:inl.find(" ")]
					glist  = glist.split("|")				       
					if glist[-1] != "" and inl.find("[") == -1:
						# PIR style header
						if glist[2] == "pir":
							glist[3] = glist[4]
							desc     = inl[inl.find(" ")+1:
										  inl.find(" - ")]
							species  = inl[inl.find(" - ")+3:]
						# GB screwup type
						elif glist[2] == "gb":
							desc     = inl[inl.find(" ")+1:]
							species  = "unknown"
						# No idea where this come form
						elif glist[2] == "sp":
							desc     = inl[inl.find(" ")+1:]   + " (" +\
									   glist[4][:glist[4].find("_")] + ")"
							species  = glist[4][glist[4].find("_")+1:]
						# all other problems
						else:
							print ("PROBLEMATIC ENTRY:")
							print (inl)		    
					else:										   
						desc   = inl[inl.find(" ")+1:inl.find("[")-1]
						species= inl[inl.find("["):]
					
					# this will be: "gi",gi,db,acc,???
					oup.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (glist[1],glist[2],
												glist[3],glist[4],desc,species))	
				# use 1st space as delimiter style
				elif style == "1stspace":
					L = inl.split(" ")
					if len(L) > 1:
						oup.write("%s\t%s\n" %
										(L[0],string.joinfields(L[1:],"\t")))
					else:
						oup.write("%s\n" % inl)
				elif style == "":
					if delim != "":
						L = string.joinfields(self.rmlb(inl).split(delim),"\t")
						oup.write("%s\n" % L)
					else:
						oup.write("%s\n" % inl)
				# all others
				else:
					L = self.rmlb(inl).split(" ")
					oup.write("%s\n" % string.joinfields(L,"\t"))
					
			inl = inp.readline()
		print ("Done!")

	def parse_ensembl_fasta(self,fasta,sp):

		print ("Parse ensembl fasta:",fasta)
		inp = open(fasta,"r")
		oup = open(fasta+".ensembl","w")
		inl = inp.readline()
		seq = ""
		c   = 0
		while inl != "":
			if inl[0] == ">":
				if seq != "":
					oup.write("%i\t%s\n" % (len(seq),seq))
					
				llist = inl[1:-1].split(" ")
				
				if len(llist) == 1:
					if llist[0].find("SINFRU") != -1:
						index = llist[0][llist[0].find("SINFRU")+6:]
						llist.append("SINFRUG"+index)					   
						llist.extend(["","","0","0",""])
					elif sp == "ce":
						try:
							int(llist[0][-1])
							llist.append(llist[0])
						except ValueError:
							llist.append(llist[0][:-1])
						llist.extend(["","","0","0",""])
					else:
						print ("Line token number mismatch, quit!")
						sys.exit(0)

				elif sp == "dm":
					tlist = []
					tlist.append(llist[1][llist[1].find(":")+1:]) # name
					tlist.append(llist[3][llist[3].find("[")+1:]) # gene
					tlist.extend(["",""])
					clist = llist[8].split(":")
					tlist.append(clist[1][1:])		    # chr
					tlist.extend(["0",""])
						
					llist = tlist
					#print llist
						
				if sp == "ce" or sp == "dm":
					id = c
				else:
					id = int(llist[0][-8:])
				
				#print llist
				oup.write("%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
					     (id,			      # id
						  llist[0],			# name
					  llist[1][llist[1].find(":")+1:], # gene
					  llist[2][llist[2].find(":")+1:], # clone
					  llist[3][llist[3].find(":")+1:], # contig
					  llist[4][llist[4].find(":")+1:], # chr
					  llist[5][llist[5].find(":")+1:], # basepair
					  llist[6][llist[6].find(":")+1:]))# status
				seq = ""
				c += 1
			else:
				seq = seq + inl[:-1]
			
			inl = inp.readline()
		print ("Done!")

	#
	# This is for getting the longest transcript or protein seq of the same
	# gene.
	#
	def ensembl_longest(self,fasta):
		# assume fasta header is in the format:
		# >pep_id pep:xxx chromosome:xxx gene:xxx transcript:xxx
		
		# put fasta into dict
		print ("Read fasta into dict:")
		inp = open(fasta)
		fd  = {}
		inl = inp.readline()
		seq = gid = pid = ""
		countG = 0
		countP = 0
		while inl != "" or seq != "":
			if inl == "" and seq != "" or inl[0] == ">":
				if seq != "":
					if gid == "" or pid == "":
						print ("Problem g or pid:",[gid,pid])
						sys.exit(0)
					if gid in fd:
						countP += 1
						if pid in fd[gid]:
							print ("This should not happen:",pid)
						else:
							fd[gid][pid] = seq
					else:
						countG += 1
						countP += 1
						fd[gid] = {pid:seq}
					# reset seq
					seq = ""
					# just put the last entry into dict, out of the loop
					if inl == "":
						break
				# set gid and pid
				L = self.rmlb(inl[1:]).split(" ")
				pid = L[0]
				for j in L[1:]:
					if j.find("gene") != -1:
						gid = j[j.find(":")+1:]
						break
				#print gid,pid
			else:
				seq += self.rmlb(inl)
			inl = inp.readline()
		
		print (" %i genes, %i peptides" % (countG,countP))
		
		# go through each gene and output the longeset one
		print ("Get the longest one...")
		oup = open(fasta+".long","w")
		for i in fd:
			longest = []
			for j in fd[i]:
				if longest == [] or longest[1] < len(fd[i][j]):
					longest = [j,len(fd[i][j])]
			oup.write(">%s\n%s\n" % (longest[0],fd[i][longest[0]]))
		
		print ("Done!")
	
	#
	# Delimiter: if "\t", should send "tab"
	#
	def simplify_desc(self,fasta,style,delim,tokens):

		print (" Process:%s" % fasta)
		print (" Output :%s.mod.fa" % fasta)
		print (" Delim  :",[delim])
		print (" Tokens :",[tokens])
		
		if delim in ["\\t","tab"]:
			delim = "\t"
		
		if tokens != "":
			tmp = tokens.split(",")
			tokens = []
			for i in tmp:
				tokens.append(int(i))
		
		print (" Style  :%s\n" % style)
		inp = open(fasta,"r")
		oup = open(fasta+".mod.fa","w")
		inl = inp.readline()
		count = 0
		while inl != "":
			inl = self.rmlb(inl)
			if inl != "" and inl[0] == ">":
				#if count % 1e3 == 0:
				#	print " %i k" % (count/1e3)
				count = count+1
				if style == "gb":
					if inl.find("|") == -1:
						print (inl)
					llist = inl.split("|")
					inl = ">" + llist[1]
				elif style == "ensembl":
					inl = inl[:inl.find(" ")]
					if inl.find(":") != -1:
						inl = ">" + inl[inl.find(":")+1:]
				elif delim != "":
					idx = inl.find(delim)
					# if delimiter does not exist, just don't change anything
					if idx == -1:
						pass
					else:
						# no token specified just take the 1st one
						if tokens == "":
							inl = inl[:idx]
						else:
							L = inl.split(delim)
							tmp = []
							for j in tokens:
								tmp.append(L[j])
							inl = "_".join(tmp)
							if inl[0] != ">":
								inl = ">"+inl
				else:
					inl = inl[:inl.find(" ")]
			
			if inl != "":		
				oup.write(inl+"\n")		       
			
			inl = inp.readline()

		oup.close()
		print (" Total %i entries" % count)
		#print "Done!"

	##
	# @param fasta   the seq file to change
	# @param prefix  species code
	##
	def add_prefix(self,fasta,prefix):
	
		print ("Add prefix...")
		inp = open(fasta,"r")
		oup = open(fasta+".mod.fa","w")
		inl = inp.readline()
		while inl != "":
			if inl[0] == ">":
				inl = ">" + prefix + "_" + inl[1:]

			if inl != "\n" and inl != "\r\n":
				if inl[-2:] == "\r\n":
					inl = inl[:-2] + "\n"
				oup.write(inl)
						
			inl = inp.readline()
			
		print ("Output file %s generated." % (fasta+".mod.fa"))
		print ("Done!")
	
	##
	# This function assumes that each line of the fasta file is ended in "\n"
	# and that no empty line is present between sequences. 
	##
	def size_filter(self,fasta,threshold,out_name):
		
		print ("Retrieve sequences with length above threshold...")
		
		if out_name == "":
			out_name = fasta + "_T%i" % threshold + ".fa"
		
		inp = open(fasta,"r")
		oup = open(out_name,"w")
		
		size = 0
		inl = inp.readline()
		llist = []
		sdict = {}
		while inl != "":		     
			if inl[0] == ">":
				# write llist if not empty and size is above threshold
				if llist != []:
					for i in llist[1:]:
						if i[-2:] == "\r\n":
							size = size + len(i) - 2
						elif i[-1] == "\n":
							size = size + len(i) - 1
						else:
							size = size + len(i)						    
					if sdict.has_key(size):
						sdict[size] = sdict[size]+1
					else:
						sdict[size] = 1
					
					if size >= threshold:
						for i in llist:
							oup.write(i)
				
				# reset values
				llist = []
				size  = 0
				
			llist.append(inl)		    
			inl = inp.readline()
		
		# deal with the last sequence
		for i in llist[1:]:
			if i[-2:] == "\r\n":
				size = size + len(i) - 2
			elif i[-1] == "\n":
				size = size + len(i) - 1
			else:
				size = size + len(i)

		if sdict.has_key(size):
			sdict[size] = sdict[size]+1
		else:
			sdict[size] = 1
		
		if size >= threshold:
			for i in llist:
				oup.write(i)
		
		print ("Fasta output %s generated" % out_name)
		
		out_name = out_name + ".stat"
		
		keys = sdict.keys()
		keys.sort()
		oup = open(out_name,"w")
		for i in keys:
			oup.write("%i\t%s\n" % (i,sdict[i]))
			
		print ("Length statistics %s generated" % out_name)
		print ("Done!")   
	
	##
	# Compare two lists and generate a output with ".comp" 
	#
	# @param files   a string of file names separated by ","
	# @param outname name for output file
	##
	def compare_lists(self,files,outname):
		
		lists = files.split(",")
		list1 = lists[0]
		list2 = lists[1]
		
		# put list1 into dict
		inp = open(list1,"r")
		dict1 = {}
		inl = inp.readline()
		while inl != "":
			if not dict1.has_key(inl[:-1].lower()):
				dict1[inl[:-1].lower()] = 0
			inl = inp.readline()
			
		# put list2 into dict
		inp = open(list2,"r")
		dict2 = {}
		inl = inp.readline()
		while inl != "":
			if not dict2.has_key(inl[:-1].lower()):
				dict2[inl[:-1].lower()] = 0  
			inl = inp.readline()
			
		
		if outname == "":
			outname = "TEMP.COMP"
		oup = open(outname,"w")
		
		# compare dict1 against dict2
		keys1 = dict1.keys()
		keys2 = dict2.keys()
		list1not2 = 0
		for i in keys1:
			if dict2.has_key(i):
				oup.write(i+"\t"+i+"\n")
				dict2[i] = 1
			else:
				list1not2 = list1not2 +1
				oup.write(i+"\t-\n")
		
		# write the rest of dict2 not compared
		list2not1 = 0	   
		for i in keys2:
			if dict2[i] == 0:
				list2not1 = list2not1 +1
				oup.write("-\t"+i+"\n") 
		
		print ("In list1 not list2:",list1not2)
		print ("   list2 not list1:",list2not1)


	##
	# The extract_cds function in ParseBlast will leave the assembly info in
	# the sequnence within the brackets. Rid of these things and replace it
	# with "X"
	##
	def cleanup(self,fasta):
		
		inp = open(fasta,"r")							   
		oup = open(fasta+"_clean.fa","w")
		
		inl = inp.readline()
		while inl != "":
			if inl.find("(") != -1:
				oup.write(inl[:inl.find("(")]+"X"+\
						  inl[inl.find(")")+1:])
			else:
				oup.write(inl)
			inl = inp.readline()


	#
	# Get rid of sequences with redundant names
	# 
	# 11/17,04 Major rewrite. Use fasta_to_dict
	#
	def del_redun_names(self,fasta):
		fdict = self.fasta_to_dict(fasta,0,verbose=1)
		fkeys = fdict.keys()
		fkeys.sort()
		oup   = open(fasta+".mod","w")
		for i in fkeys:
			oup.write(">%s\n%s\n" % (i,fdict[i]))
		
		print ("Done!")
	
	#
	# Delete sequences that are completely the same
	# @delete    Rid of the redun ones (1, default) or just generate a list (0)
	# @wsp       Only consider within species (0) or don't care (0, default)
	# @sp_prefix Length of sp header to consider
	#
	def del_redun_seqs(self,fasta,delete=1,wsp=0,sp_prefix=0):
		F  = self.fasta_to_dict(fasta,0,verbose=1)
		fk = F.keys()
		fk.sort()
		R  = {} # {redundant_with_name:[redundant]}
		S  = {} # {non_redun_name:non_redun_seq}
		for i in range(len(fk)):
			si = F[fk[i]]
			if wsp:
				spi = fk[i][:sp_prefix]  # species for seq i
			if si not in R:
				S[fk[i]] = F[fk[i]]
			for j in range(i+1,len(fk)):
				sj = F[fk[j]]
				if wsp:
					spj = fk[j][:sp_prefix]# species for seq j
				# same sequence
				if si == sj:
					q = 1
					# consider only within sp and not from the same sp
					if wsp and spi != spj:
						q = 0
					if q:
						if fk[i] not in R:
							R[fk[i]] = [fk[j]]
						else:
							R[fk[i]].append(fk[j])

		print ("Redundant:")
		for i in R:
			print (i,len(R[i]),R[i])
		
		if delete:
			oup = open(fasta+".mRedunSeq","w")					
			for i in S:
				oup.write(">%s\n%s\n" % (i,S[i]))
	
		oup = open(fasta+".RedunSeq_names","w")
		for i in R:
			oup.write("%s\n" % i)

		oup.close()
		print ("%i seq redun" % len(R.keys()))
		print ("Done!")
	
	def check_redun(self,dir):		
		f = os.listdir(dir)
		for i in f:
			try:	
				print (i)
				self.fasta_to_dict(i,0,verbose=1)
				print ("")
			except IndexError:
				continue
		
		print ("Done!")


	##
	# This is written to get rid of sequences within a particular fasta file
	#
	# @param fasta  fasta file with entries to be deleted
	# @param dlist  file with a list of names to be deleted, or a string with
	#               species header. ALWAYS look for 2 charcters
	##
	def delete(self,fasta,dlist):
		
		dflag = 0
		try:
			print ("Read delete list:")
			ddict = f_util.file_to_dict(dlist,0)
			print ("To be deleted:",len(ddict.keys()))
		except IOError:
			print (" not a file.")
			dlist = dlist.split(",")
			dflag = 1
		
		
		inp = open(fasta,"r")
		oup = open(fasta+".mod","w")
		inl = inp.readline()
		write  = 1
		countD = 0
		countW = 0
		countT = 0
		while inl != "":
			inl = self.rmlb(inl)
			if inl == "":
				inl = inp.readline()
				continue
			
			if inl[0] == ">":
				countT += 1					
				if dflag:
					if inl[1:3] in dlist:
						write = 0
						countD += 1
					else:
						write = 1
						countW += 1
				else:
					if ddict.has_key(inl[1:]):
						ddict[inl[1:]] = -1
						write = 0
						countD += 1
					else:
						write = 1
						countW += 1
			if write == 1:
				oup.write(inl+"\n")
			
			inl = inp.readline()
		
		if dflag == 0:
			oup = open(fasta+".log","w")
			for i in ddict:
				if ddict[i] != -1:
					oup.write("%s\n" % i)
		
		print ("Total %i, write %i, deleted %i" % (countT,countW,countD))
		print ("Deletion done!")
		
	#
	# Break fasta file into several different files.
	#  04/13,10 Last file too big. Fix it.
	# 
	# @param by     divide by the int passed
	# @param setdir make dir for each file
	#	
	def divide(self,fasta,by,setdir=0,verbose=0,newline=0):
		
		if by == 1:
			print ("What's the point divided by 1?? Quit!")
			sys.exit(0)
		
		fdict = self.fasta_to_dict(fasta,0,verbose,newline)
		fkeys = sorted(fdict.keys())
		#fkeys.sort()
		step  = int(len(fkeys)/by+1)
		
		c = 1 # increment, file number
		print ("Number of seq   :",len(fkeys))
		print ("Avg size divided:",step)
		for i in range(0,len(fkeys),step):
			oup  = open(fasta+"_%i" % c,"w")
			keys = fkeys[i:i+step]
			#print len(keys), keys[0],keys[-1]
			
			for j in keys:
				oup.write(">%s\n%s\n" % (j,fdict[j]))
			oup.close()
			
			if setdir:
				os.system("mkdir %i" % i)
				os.system("mv %s_%i %i" % (fasta,i,i))			
			c += 1
		
		"""
		C = 0
		for i in range(by):
			oup  = open(fasta+"_%i" % (i+1),"w")
			if C == by -1:
				keys = fkeys[step*C:]	
			else:
				keys = fkeys[step*C:step*(C+1)]
			print len(keys), keys[0],keys[-1]
			
			for j in keys:
				oup.write(">%s\n%s\n" % (j,fdict[j]))
			
			if setdir:
				os.system("mkdir %i" % i)
				os.system("mv %s_%i %i" % (fasta,i,i))
			C = C+1
		"""	
		print ("Done!")
		
	def divide1seq(self,fasta,size):
		inp = open(fasta)
		sid = inp.readline()
		if sid[0] != ">":
			print ("Fasta file?", fasta)
			print ("Quit!")
			sys.exit(0)
		
		sid = sid[1:].strip()
		seq = ""
		inl = inp.readline()
		while inl != "":
			seq += inl.strip()
			inl = inp.readline()
		
		oup = open(fasta+".divided_%i" % size,"w")
		seqL= len(seq)
		frag= seqL/size + 1
		print ("Length:",seqL)
		for i in range(0,seqL,size):
			print ("",i+1,i+size)
			oup.write(">%s|%i-%i\n%s\n" % (sid,i+1,i+size,seq[i:i+size]))
		
		print ("Done!")
		

	#
	# Concatenate sequences in a file together and record the coordinates
	#
	def concat(self, fasta):
		
		fdict = self.fasta_to_dict(fasta)
		oup1  = open(fasta+".concat","w")
		oup1.write(">%s.concat\n" % fasta)
		oup2  = open(fasta+".concat.coord","w")
		
		c = 0
		for i in fdict:
			seq = self.rmlb(fdict[i][1])
			
			oup1.write("%s" % seq)
			oup2.write("%s\t%i\t%i\n" % (i,c+1,c+len(seq)))
			c = c + len(seq)
		
		print ("Done!")
	
	#
	# Extension of concat. Use in the following situation, for example:
	#
	# Several rice contig sequences are concatenated together. Then genes are
	# predicted from such concatenated sequence file. Now we want to map the
	# genes back to the contigs.
	#
	# @param coord1 the coord output file of concat. L should be smaller than R.
	# @param coord2 a file with [seq][L][R]. The coordinates should be within
	#               the range of coord1. L does not have to be smaller than R.
	#
	def concat_locate(self, coord1, coord2):
		
		c1  = f_util.file_to_dict(coord1,3)
		c2  = f_util.file_to_dict(coord2,3)
		oup = open(coord2+".mapped","w")
		for i in c2:
			if int(c2[i][0]) < int(c2[i][1]):
				gL = int(c2[i][0])
				gR = int(c2[i][1])
			else:
				gL = int(c2[i][1])
				gR = int(c2[i][0])
			
			hit = "not_mapped"
			for j in c1:
				#print int(c1[j][0]),int(c1[j][1]),"--->",gL,gR
				if int(c1[j][0]) < gL and int(c1[j][1]) > gR:
					hit = j
					#print "HERE"
					break
			oup.write("%s\t%s\n" % (i,hit))
			
		print ("Done!")
	
	#
	# Store each fasta squence in a file into its own file
	#
	def indiv(self,fasta,odir):
		
		if odir != "":
			os.system("mkdir %s" % odir)
		
		inp = open(fasta)
		inl = inp.readline()
		print ("Write seq...")
		c = 0
		while inl != "":
			if inl[0] == ">":
				if c % 1e3 == 0:
					print (" %i k" % (c/1e3))
				c += 1
				seq = self.rmlb(inl[1:])
				seq = string.joinfields(seq.split(" "),"_")
				oup = open("%s/%s.fa" % (odir,seq),"w")
				oup.write(">%s\n" % seq)
			else:
				oup.write(inl)
			inl = inp.readline()
		
		print ("Done!")
	
	def count(self,fasta):
		inp = open(fasta)
		inl = inp.readline()
		print ("Count seq...")
		c = 0
		while inl != "":
			if inl[0] == ">":
				c += 1
			inl = inp.readline()
		print ("Total %i sequences" % c)
		print ("Done!")	

	#
	# @param D the delimiter for alt spliced designation, only consider the last
	#
	def get_longest(self,fasta,D):
		
		fa = self.fasta_to_dict(fasta)
		fk = fa.keys()
		gn = {} # gn = {gene_name:{ver:length}}
		countT = 0
		countG = 0
		oup = open("log_get_longest","w")
		oup.write("Gene name format problem:\n")
		for i in fk:
			countT += 1
			g = i.split(D)		# gene
			#print g
			# No deliminter
			if len(g) < 2:
				oup.write(g[0]+"\n")
				N = g[0]
				A = ""
			else:
				# Only consider the last token for alt #
				N = D.join(g[:-1])	# name without alt #
				A = g[-1]			# alt #
				
			if N not in gn:
				countG += 1
				gn[N] = {A:len(fa[i])}
			else:
				gn[N][A] = len(fa[i])
				
		
		oup1 = open(fasta+".longest","w")
		oup2 = open(fasta+".longest.log","w")
		for i in gn:
			# compare length of a gene_name
			maxLen = 0
			maxVer = ""
			for j in gn[i]:
				if gn[i][j] > maxLen:
					maxLen = gn[i][j]
					maxVer  = j
			if maxVer != "":
				maxID = "%s%s%s" % (i,D,maxVer)
			else:
				maxID = i
			oup1.write(">%s\n%s\n" % (maxID,fa[maxID]))
			for j in gn[i]:
				flag = "-"
				if j == maxVer:
					flag = "long"
				ID = "%s%s%s" % (i,D,j)
				oup2.write("%s\t%s\t%i\t%s\n" % (i,j,gn[i][j],flag))
				
		oup.close()
		oup1.close()
		oup2.close()
		
		print (" Total %i pep belong to %i genes" % (countT,countG))
	
	# Calculat the GC content
	def get_gc(self,fasta):
		fa = self.fasta_to_dict(fasta)
		GC = []
		for i in fa:
			s = fa[i]
			g = len(s.split("G"))-1
			c = len(s.split("C"))-1
			GC.append(float(g+c)/len(s))
			print (i,float(g+c)/len(s))
		
		print ("GC content:",sum(GC)/len(GC))
	
	def get_names(self,fasta,call=0):
		fa = self.fasta_to_dict(fasta)
		fk = fa.keys()
		fk.sort()
		
		if call:
			return fk
		else:
			fk = "\n".join(fk)
			oup = open(fasta+".names","w")
			oup.write(fk)
			print ("Done!")
	
	#
	# Format fasta file into certain width and delete stop char if asked
	#
	def format(self,fasta,linew,ridstop):
		
		fdict = self.fasta_to_dict(fasta)
		oup = open("%s.formatted" % fasta,"w")
		for i in fdict:
			oup.write(">%s\n" % i)
			seq = fdict[i]
			if ridstop and seq[-1] == "*":
				seq = seq[:-1]
			
			for i in range(0,len(seq),linew):
				oup.write("%s\n" % seq[i:i+linew])
		
		print ("Done!")
			
	def rmlb(self,astr):
		if astr[-2:] == "\r\n":
			astr = astr[:-2]
		elif astr[-1] == "\n":
			astr = astr[:-1]
		return astr
		

	def help(self):
		print (" -f function")
		print ("    get_sequences - get specified sequences. REQUIRES: -fasta")
		print ("       -name, Optional: -type, -match")
		print ("    getseq2 - simpler function. REQUIRES: fasta,name. OPT:tokens")
		print ("    get_group_seq - get sequences for groups based on the group")
		print ("       specification passed. REQUIRES: -fasta, -group. Optional:")
		print ("       create_d, includ_o")
		print ("    get_stretch2 - get segment of a SINGLE sequence. Best for")
		print ("       chromosome seq. REQUIRES: fasta, coords, OPT: isfile")
		print ("    get_stretch3 - Take single seq fasta REQUIRES: fasta, coords")
		print ("    get_stretch4 - Take multiple seq fasta, REQ: fasta,coords") # use get_stretch4 where you have Chr\tseq1,seq2\tname
		print ("    get_gc - calculate GC content")
		print ("    gff_to_coord - convert GFF to coord file for get strech,")
		print ("       NEED: gff")
		print ("    gff_promoter_to_coord - convert GFF to promoter coord file")
		print ("       NEED: gff")
		print ("    gff_prom_to_coord2 - convert GFF to promoter coord file where")
		print ("       promoter is not specified (get 1000bp upstream) - NEED: gff")
		print ("    gff_prom_to_coord_5utr - convert GFF to promoter coord file where")
		print ("      promoter is not specified (get 1000bp upstream and first 200bp of")
		print ("      5'UTR) - NEED: gff")
		print ("    gff_cds_to_coord - convert GFF to cds coord file - NEED: gff")
		print ("    fasta_to_oneline - convert fasta file to one line")
		print ("       format. Requires: -fasta. OPT: d")
		print ("    fasta_to_phylip - NEED: fasta, OPT: nlen")
		print ("    fasta_to_stockholm - NEED: fasta")
		print ("    oneline_to_fasta - convert the other way around")
		print ("	      REQUIRES: -oneline")
		print ("    parse_desc - get the description out of fasta file")
		print ("       NEED: fasta, OPT: style, D")
		print ("    simplify - simplify the header by taking the id imme-")
		print ("       diately after '>'.")
		print ("       Requires: -fasta, Optional:-style, -D, tokens")
		print ("    prefix - add whatever prefix to all sequences in")
		print ("       a fasta file. Requires -fasta, -prefix")
		print ("    size_filter - get sequences longer than the treshold")
		print ("       passed. Requires: -fasta, -T. Optional: out_name")
		print ("    compare_lists - compare two lists, requires: -files,")
		print ("       optional: out")
		print ("    cleanup - get rid of bracketed characters and replace")
		print ("       with an '-'. Requires: -fasta")
		print ("    get_sizes - get the sizes of sequences within a fasta")
		print ("    divide - divide fasta file into parts. NEED: fasta,by,setdir")
		print ("       file. REQUIRES: -fasta. OPTIONAL: x, newline")
		print ("    divide1seq - divide 1 sequence in a fasta. NEED: fasta, size")
		print ("    index_names - convert names into indices. Generate a")
		print ("       fasta file with seq_id converted to index, and a")
		print ("       '.name' file with both seq_id and indices")
		print ("       REQUIRES: -fasta")
		print ("    index_pairs - convert names from a pair file into indices")
		print ("       NEED: fasta, pairs")
		print ("    change_names - change names in the Fasta header.")
		print ("       REQUIRES: -fasta, -name, -desc OPTIONAL:-delim")
		print ("    rename - Change id within the fasta file")
		print ("       REQUIRES: -fasta, -name")
		print ("    rename_all - change id of all .fa or .fasta file in a dir")
		print ("       NEED: -targetd, -name")
		print ("    del_redun_names - delete seq with redundant names. NEED: -fasta")
		print ("    del_redun_seqs - delete identical sequences and generate a list,")
		print ("       NEED: -fasta, OPT: delete")
		print ("    check_redun - check fasta in a folder. NEED: targetd")
		print ("    get_sp - get species for a particular species")
		print ("       REQUIRES: -fasta, -sp")
		print ("    count_sp - count # of seq. NEED: fasta, sp")
		print ("    parse_ensembl_fasta - parse desc, get len, one line")
		print ("       seqeunce. REQUIRES: -fasta (in ensembl format)")
		print ("    covert_header - coord from id<space>L-R to id_((L-1)*3+1)_R,")
		print ("       NEED: fasta")
		print ("    delete - rid of seqeunces, REQUIRES: -fasta,-dlist")
		print ("    concat - concatenate sequences in a file, REQUIRES: -fasta")
		print ("    locate - a companion function of concat. REQUIRES: c1, c2")
		print ("    indiv  - Store each fasta squence in a file into its own")
		print ("        REQUIRES: fasta, OPT: odir")
		print ("    mask   - masking the areas based on the passed coord.")
		print ("        REQUIRES: fasta, coords")
		print ("    count  - count the number of sequences, REQUIRES: fasta")
		print ("    get_longest - get the longest among alternative predictions")
		print ("        NEED: fasta, D")
		print ("    format - format fasta line into certain width and rid of")
		print ("        the stop character if choose to. NEED: fasta, OPT: ")
		print ("        linew, ridstop")
		print (" -fasta : location and name of fasta file for getting or con-")
		print ("	  verting sequences.")
		print (" -group : group specification [group_name][seq_id]")
		print (" -oneline: sequence in the format [name][seq]")
		print (" -files : names separated by ','" )
		print (" -name  : file with lines of gene names. For change_names(), the")
		print ("	  file is in the format [new_name][old_name]")
		print (" -newline: add new line char")
		print (" -nlen  : name length threshold, for fasta_to_phylip")
		print (" -seqid : (1) use 4th col in coord file as seqid or not (0, default)")
		print (" -style : The fasta description file format:")
		print ("	    GenBank - gb")
		print (" -size  : size of sequence blocks to divide into")
		print (" -prefix: the prefix to be added to the sequence descriptors" )
		print (" -out   : output file name")
		print (" -T     : threshold length of sequences included in the output")
		print (" -start : the starting index, default = 1")
		print (" -type  : peptide [0,default], nucleotide [1]")
		print (" -D     : delimiter")
		print (" -species: two CAPITAL characters species abbreviation")
		print (" -dlist : a list of names to be deleted from fasta file")
		print (" -delim : delimiter between id and other info")
		print (" -match : 0: exacth match to fasta descriptor (default) or 1:")
		print ("	  exact match to the passed name right after '>'")
		print (" -desc  : whether the name contain description [1, default] or")
		print ("          not [0]")
		print (" -coords: the coordinates of the segments to get, the number of")
		print ("	  coordinates should be even and none of them should be")
		print ("	  smaller than 1. They are separated by ','. Alternatively,")
		print ("   For get_stretch2: file contains two coord separated by ',' ")
		print ("       get_stretch3: whatever...[L][R]")
		print ("       get_stretch4: two column file with seq name and coord, coord")
		print ("          in the format L1,R1,... or two coords separated by ','")
		print ("          passed through command line")
		print (" -delete: Delete redun seq (1, default) or not (0)")
		print (" -wcoord: output sequence name with coords")
		print (" -create_d: create dir for each sequence groups, default no [0]")
		print (" -includ_o: include outgroup sequence [1], or not [0, default]")
		print ("          For this to work, the fasta file passed should have a")
		print ("          sequence called OUT.")
		print (" -ignore: For rename, ignore the string after the specified char")
		print ("          during name matching. And add the string back later")
		print (" -format: 0 (default) for pep or nt coords for pep or nt seq file")
		print ("          respectively. 1 for pep coords to get nt seq.")
		print (" -verbose: count line for fasta_to_dict [1], or not [0, default")
		print (" -by    : divide fasta file by the passed integer")
		print (" -add   : add domain info in get_strech. The last token in the")
		print ("          coord file should be domain name. Default 0")
		print (" -x     : count masked region size")
		print (" -wsp   : consider within sp (1) or don't care (0, default)")
		print (" -tokens: indice seperated by ',' defining the tokens to be")
		print ("          taken as names in the passed name file. DEFAULT ''")
		print (" -setdir: make dir for each divided file [1], default no [0]")
		print (" -sp_prefix lenght of species prefix")
		print (" -targetd : target dir")
		print (" -d     : delimiter for fasta_to_oneline, defailt ''")
		print (" -linew : format line width, default 80")
		print (" -ridstop: rid of stop character, default no [0]")
		print (" -odir  : output directory")
		print ("")
		sys.exit(0)

		
#
# Function call
#

f_util = FileUtility.file_util()

if __name__ == '__main__':
	fasta = name_file = function = style = prefix = out_name = files = \
			oneline   = delimiter= species = dlist = delim = coords = group = \
			ignore    = by = c1 = c2 = tokens = targetD = d = odir = size = \
			gff = pairs = "" 
	threshold = type = match = create_d = includ_o = format = verbose = \
			wcoord = segment = add =  call = x = setdir = ridstop = newline = nlen = \
			sp_prefix = seqid = 0
	start = desc = isfile = delete = xsp = 1
	linew     = 80
	manager   = fasta_manager()
	trans     = Translation.translate()

	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-fasta":
			fasta      = sys.argv[i+1]
		elif sys.argv[i] == "-oneline":
			oneline    = sys.argv[i+1]
		elif sys.argv[i] == "-name":	      
			name_file  = sys.argv[i+1]
		elif sys.argv[i] == "-nlen":	      
			nlen       = int(sys.argv[i+1])
		elif sys.argv[i] == "-f":	 
			function   = sys.argv[i+1]
		elif sys.argv[i] == "-seqid":	     
			seqid      = int(sys.argv[i+1])
		elif sys.argv[i] == "-style":	     
			style      = sys.argv[i+1]
		elif sys.argv[i] == "-prefix":	    
			prefix     = sys.argv[i+1]
		elif sys.argv[i] == "-out":	       
			out_name   = sys.argv[i+1]
		elif sys.argv[i] == "-files":	     
			files      = sys.argv[i+1]	      
		elif sys.argv[i] == "-D":	 
			delimiter  = sys.argv[i+1]	      
		elif sys.argv[i] == "-T":	 
			threshold  = int(sys.argv[i+1])
		elif sys.argv[i] == "-start":	     
			start      = int(sys.argv[i+1])
		elif sys.argv[i] == "-segment":	   
			segment    = int(sys.argv[i+1])
		elif sys.argv[i] == "-type":	      
			type       = int(sys.argv[i+1])
		elif sys.argv[i] == "-match":	     
			match      = int(sys.argv[i+1])
		elif sys.argv[i] == "-species":	   
			species    = sys.argv[i+1]
		elif sys.argv[i] == "-dlist":	     
			dlist      = sys.argv[i+1]
		elif sys.argv[i] == "-delim":	     
			delim      = sys.argv[i+1]
		elif sys.argv[i] == "-delete":	      
			delete     = int(sys.argv[i+1])
		elif sys.argv[i] == "-desc":	      
			desc       = int(sys.argv[i+1])
		elif sys.argv[i] == "-coords":	    
			coords     = sys.argv[i+1]
		elif sys.argv[i] == "-group":	    
			group      = sys.argv[i+1]
		elif sys.argv[i] == "-create_d":	     
			create_d   = int(sys.argv[i+1])
		elif sys.argv[i] == "-includ_o":	     
			includ_o   = int(sys.argv[i+1])
		elif sys.argv[i] == "-ignore":	    
			ignore     = sys.argv[i+1]
		elif sys.argv[i] == "-format":	    
			format     = int(sys.argv[i+1])			
		elif sys.argv[i] == "-verbose":	    
			verbose    = int(sys.argv[i+1])			
		elif sys.argv[i] == "-by":	    
			by         = int(sys.argv[i+1])			
		elif sys.argv[i] == "-c1":	    
			c1         = sys.argv[i+1]
		elif sys.argv[i] == "-c2":	    
			c2         = sys.argv[i+1]
		elif sys.argv[i] == "-wcoord":	    
			wcoord     = int(sys.argv[i+1])
		elif sys.argv[i] == "-add":	    
			add        = int(sys.argv[i+1])
		elif sys.argv[i] == "-call":	    
			call        = int(sys.argv[i+1])
		elif sys.argv[i] == "-isfile":	    
			isfile      = int(sys.argv[i+1])
		elif sys.argv[i] == "-x":	    
			x           = int(sys.argv[i+1])
		elif sys.argv[i] == "-tokens":	    
			tokens      = sys.argv[i+1]
		elif sys.argv[i] == "-targetd":
			targetD     = sys.argv[i+1]
		elif sys.argv[i] == "-d":
			d     = sys.argv[i+1]
		elif sys.argv[i] == "-linew":
			linew = int(sys.argv[i+1])
		elif sys.argv[i] == "-ridstop":
			ridstop = int(sys.argv[i+1])
		elif sys.argv[i] == "-odir":
			odir = sys.argv[i+1]
		elif sys.argv[i] == "-size":
			size = int(sys.argv[i+1])
		elif sys.argv[i] == "-newline":
			newline = int(sys.argv[i+1])
		elif sys.argv[i] == "-gff":
			gff = sys.argv[i+1]
		elif sys.argv[i] == "-pairs":
			pairs = sys.argv[i+1]
		elif sys.argv[i] == "-wsp":
			wsp = int(sys.argv[i+1])
		elif sys.argv[i] == "-sp_prefix":
			sp_prefix = int(sys.argv[i+1])
		else:
			print ("Unknown parameter:",sys.argv[i])
			print ("Quit!")
			sys.exit(0)

	if function == "get_sequences":
		if fasta == "" or name_file == "":
			print ("\nNeed fasta file and name file\n")
			manager.help()
		manager.get_sequences(fasta,name_file,segment,type,match)	
	elif function == "getseq2":
		if fasta == "" or name_file == "":
			print ("\nNeed fasta and name.\n")
			manager.help()
		manager.getseq2(fasta,name_file,tokens)
	elif function == "get_group_seq":
		if fasta == "" or group == "":
			print ("\nNeed fasta file and group specification.\n")
			manager.help()
		manager.get_group_seq(fasta,group,create_d,includ_o)
	elif function == "get_stretch":
		if fasta == "" or coords == "":
			print ("\nNeed fasta file and coordinates\n")
			manager.help()
		manager.get_stretch(fasta,coords,format,wcoord,verbose,call,add)
	elif function == "get_stretch2":
		if fasta == "" or coords == "":
			print ("\nNeed fasta file and coordinates\n")
			manager.help()
		manager.get_stretch2(fasta,coords,call,isfile)
	elif function == "get_stretch3":
		if fasta == "" or coords == "":
			print ("\nNeed fasta file and coordinates\n")
			manager.help()
		manager.get_stretch3(fasta,coords)
	elif function == "get_stretch4":
		if fasta == "" or coords == "":
			print ("\nNeed fasta file and coordinates\n")
			manager.help()
		manager.get_stretch4(fasta,coords,seqid)
	elif function == "get_names":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.get_names(fasta)
	elif function == "get_gc":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.get_gc(fasta)
	elif function == "gff_to_coord":
		if gff == "":
			print ("\nNeed gff file\n")
			manager.help()
		manager.gff_to_coord(gff)
	elif function == "gff_promoter_to_coord":
	        if gff == "":
	            print ("\nNeed gff file\n")
	            manager.help()
	        manager.gff_promoter_to_coord(gff)
	elif function == "gff_prom_to_coord2":
	        if gff == "":
	           print ("\nNeed gff file\n")
	           manager.help()
	        manager.gff_prom_to_coord2(gff) 
	elif function == "gff_prom_to_coord_5utr":
	        if gff == "":
	           print ("\nNeed gff file\n")
	           manager.help()
	        manager.gff_prom_to_coord_5utr(gff)
	elif function == "gff_cds_to_coord":
			if gff == "":
				print ("\nNeed gff file\n")
				manager.help()
			manager.gff_cds_to_coord(gff)
	elif function == "fasta_to_oneline":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.fasta_to_oneline(fasta,d)
	elif function == "fasta_to_phylip":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.fasta_to_phylip(fasta,nlen)
	elif function == "fasta_to_stockholm":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.fasta_to_stockholm(fasta)
	elif function == "oneline_to_fasta":
		if oneline == "":
			print ("\nNeed one line file\n")
			manager.help()
		manager.oneline_to_fasta(oneline)
	elif function == "simplify":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.simplify_desc(fasta,style,delimiter,tokens)
	elif function == "prefix":
		if fasta == "" or prefix == "":
			print ("\nNeed fasta file and prefix\n")
			manager.help()
		manager.add_prefix(fasta,prefix)
	elif function == "size_filter":
		if fasta == "" and threshold == 0:
			print ("\nNeed fasta files and length threshold\n")
			manager.help()
		manager.size_filter(fasta,threshold,out_name)
	elif function == "parse_desc":
		if fasta == "":
			print ("\nNeed fasta file and description style info\n")
			manager.help()
		manager.parse_desc(fasta,style,delimiter)
	elif function == "compare_lists":
		if files == "" or len(files.split(",")) != 2:
			print ("\nNeed files and it should be exactly 2 names\n")
			manager.help()
		manager.compare_lists(files,out_name)
	elif function == "cleanup":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.cleanup(fasta)  
	elif function == "get_sizes":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.get_sizes(fasta,0,x)
	elif function == "index_names":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.index_names(fasta,start)
	elif function == "index_pairs":
		if fasta == "" or pairs == "":
			print ("\nNeed fasta and pair files\n")
			manager.help()
		manager.index_pairs(fasta,pairs)
	elif function == "change_names":
		if fasta == "" or name_file == "" or desc == -1:
			print ("\nNeed fasta file and name file. Also need to specify desc\n")
			manager.help()
		manager.change_names(fasta,name_file,desc,delim)
	elif function == "rename":
		if fasta == "" or name_file == "":
			print ("\nNeed fasta and name files\n")
			manager.help()
		manager.rename2(fasta,name_file,ignore) 
	elif function == "rename_all":
		if "" in [targetD,name_file]:
			print ("\nNeed target dir and name files\n")
			manager.help()
		manager.rename_all(targetD,name_file)
	elif function == "del_redun_names":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.del_redun_names(fasta)
	elif function == "del_redun_seqs":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.del_redun_seqs(fasta,delete,wsp,sp_prefix)
	elif function == "check_redun":
		if targetD == "":
			print ("\nNeed dir\n")
			manager.help()
		manager.check_redun(targetD)
	elif function == "convert_header":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.convert_header(fasta)
	elif function == "get_sp":
		if fasta == "" or species == "":
			print ("\nNeed fasta file and a species name specified\n")
			manager.help()
		manager.get_sp(fasta,species)
	elif function == "count_sp":
		if fasta == "" or species == "":
			print ("\nNeed fasta file and a species name specified\n")
			manager.help()
		manager.count_sp(fasta,species)
	elif function == "parse_ensembl_fasta":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.parse_ensembl_fasta(fasta,species)
	elif function == "delete":
		if fasta == "" or dlist == "":
			print ("\nNeed fasta file and list of names\n")
			manager.help()
		manager.delete(fasta,dlist)
	elif function == "divide":
		if fasta == "" or by == "":
			print ("\nNeed fasta file and the number of divisions\n")
			manager.help()
		manager.divide(fasta,by,setdir,verbose,newline)
	elif function == "divide1seq":
		if fasta == "" or size == "":
			print ("\nNeed fasta file and the size of sequenec blocks\n")
			manager.help()
		manager.divide1seq(fasta,size)
	elif function == "concat":
		if fasta == "":
			print ("\nNeed fasta file\n")
			manager.help()
		manager.concat(fasta)
	elif function == "locate":
		if c1 == "" or c2 == "":
			print ("\nNeed coordinate files\n")
			manager.help()
		manager.concat_locate(c1,c2)
	elif function == "indiv":
		if fasta == "":
			print ("\nNeed fasta\n")
			manager.help()
		manager.indiv(fasta,odir)
	elif function == "mask":
		if fasta == "" or coords == "":
			print ("\nNeed fasta and coords\n")
			manager.help()
		manager.mask(fasta,coords)
	elif function == "fasta_to_dict":
		manager.fasta_to_dict(fasta)
	elif function == "count":
		if fasta == "":
			print ("\nNeed fasta\n")
			manager.help()
		manager.count(fasta)
	elif function == "ensembl_longest":
		if fasta == "":
			print ("\nNeed fasta\n")
			manager.help()
		manager.ensembl_longest(fasta)
	elif function == "get_longest":
		if "" in [fasta,delimiter]:
			print ("\nNeed fasta and deliminter\n")
			manager.help()
		manager.get_longest(fasta,delimiter)
	elif function == "format":
		if "" in [fasta]:
			print ("\nNeed fasta\n")
			manager.help()
		manager.format(fasta,linew,ridstop)
	else:
		print ("\nUnknown function...\n")
		manager.help()

"""
		print "    get_stretch - get segments of all sequences in fasta"
		print "       file, this is meant for alignment output. Need:"
		print "       -fasta, -coords. OPTIONAL: -format, -wcoord, -verbose"

	##
	# If multiple coords are passed, only the last one of the coords for each
	# seqeunce will be stored. Implement multiple get later...
	#
	# @param coords  this can be a string with coords separated by "," or a file
	#                with coordinates specified.
	# @param format  0 (default). 1 to get nt sequences based on peptide coords.
	#                The sequence file passed HAS TO BE nt and synchronized with
	#                the peptide sequences where coords are derived from.
	# @param verbose for fasta_to_dict
	# @param wc      ouput sequence name with coordinate or not
	# @param call    non-command line call [1]. Return a dict instead of
	#                file output.
	# @param add     add domain name to the sequence name. The coords passed
	#                have to be [gene][L][R][dom]
	##
	def get_stretch(self,fasta,coords,format=0,wc=0,verbose=0,call=0,add=0):
		
		print "Fasta  :",fasta
		print "Coord  :",coords
		print "Format :",format
		print "W/coord:",wc
		print "verbose:",verbose
		print "Call   :",call
		print "Add_dom:",add
		
		if format:
			print "CAUTION: you try to get nt sequences using pep coords"
		
		print "Read sequence into dict"
		fdict = self.fasta_to_dict(fasta,0,verbose)
		fkeys = fdict.keys()
		print "%i sequences..." % len(fkeys)
		
		#print fdict	
		cdict = {}
		only1 = 0  # all sequences have the same coords
		
		print "Read coords..."
		# internal call, dict will be passed
		if call == 1:
			cdict = coords
		else:
			# check if coords is a file
			try:
				inp = open(coords,"r")
				inl = inp.readline()
				while inl != "":
					L = self.rmlb(inl).split("\t")
					if cdict.has_key(L[0]):
						cdict[L[0]].append(L[1:])					
					else:
						cdict[L[0]] = [L[1:]]
					inl = inp.readline()
			except IOError:
				print "Only 1 pair of coords..."
				only1 = 1
				coords = coords.split(",")
				for i in fdict:
					cdict[i] = 0
				
		# for each sequence
		oup = open(fasta+".seg.fa","w")
		count_g = 0
		count_t = 0
		count_a = 0
		ckeys = cdict.keys()
		ckeys.sort()
		for i in ckeys:
			#print i
				print "E1: seq absent,",i
				count_a += 1
				continue
			seq = fdict[i]
			slen = len(seq)
			out3prime = 0
			
			if only1 and len(coords) != 2:
				# take multiple stretches
				c = 0
				seg = ""
				while c < len(coords):
					try:
						if int(coords[c]) < 1 or int(coords[c+1]) < 1:
							print " coords smaller than 1, quit!\n"
							sys.exit(0)
					except ValueError:
						print "valueerr:",coords
						sys.exit(0)
					
					coords[c]   = int(coords[c])
					coords[c+1] = int(coords[c+1])
					if format:
						coords[c]   = coords[c]*3-2
						coords[c+1] = coords[c+1]*3

					seg = seg + seq[coords[c]-1:coords[c+1]]			    
					if coords[c+1] > slen:
						out3prime = 1			
					c += 2
					
				if wc:
					oup.write(">%s_%s\n%s\n"%(i,string.joinfields(coords),seg))
				else:
					oup.write(">%s\n%s\n" % (i,seg))
				count_t += 1
			else:
				if only1:
					clist = [coords]
				else:
					clist = cdict[i]
				count_c = 1
				
				for j in clist:
					seg = ""
					
					# check the first pair, if it is in reverse ori and there
					# are more than one pair of coord, will stop right here.
					if len(j) > 2 and int(j[0]) > int(j[1]):
						print "Coord format problem... Quit!"
						sys.exit(0)
					
					# iterate coordinate sets, some may be multiple
					for k in range(0,len(j)-1,2):
						cL = int(j[k])
						cR = int(j[k+1])
						if int(j[k]) > int(j[k+1]):
							cR = int(j[k])
							cL = int(j[k+1])
						if format:
							cL = cL*3-2
							cR = cR*3						
						seg += seq[cL-1:cR]
											
					# add domain name, dom should be last token
					if add:
						idx = "%s|%s" % (i,j[-1]) 
					elif len(clist) == 1:
						idx = i
					else:
						idx = "%s_seg%i" % (i,count_c)
					if wc:
						idx += "_%s" % (string.joinfields(j,"|"))

					# now, the coord may be in reverse. only check the first
					# pair to see if this is the case. if so, the reverse
					# complement is generated. Obviously won't apply to peptide
					# sequences
					if int(j[0]) > int(j[1]):
						seg = trans.rc(seg)
					
					oup.write(">%s\n%s\n" % (idx,seg))
					count_c += 1
					
					if int(j[1]) > slen:
						out3prime = 1
					count_t += 1
			
			count_g += 1
			if out3prime:
				pass
				#print "%s:out of 3' end, take as much as it can" % i
		
		print "Total %i stretches in %i genes, %i id not in seq file" % \
				(count_t,count_g,count_a)
		print "Done!"					
"""
