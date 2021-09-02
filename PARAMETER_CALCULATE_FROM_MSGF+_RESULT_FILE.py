#!/usr/bin/python
#NKS Program for read csv file and generate more parameter for machine learning
#3/4/2019
import csv
import re


with open("Mix3_LTQFT_with_oxidation_and_mc2.csv", 'r', newline='') as seq, open('result.csv','w',newline='') as csvFile, open ("Mix3_LTQFT.mgf",'r') as fh:
	csv_reader = csv.reader(seq, delimiter=',')
	csv_reader1 = next(csv_reader)

	x = 0
	d = dict()
	line =fh.readline()
	while line:
		line = line.rstrip('\n')
		#print('\n')
		if 'BEGIN IONS' in line:
			#d[key] = key_value
			key = x
			d[key] = ''
			x +=1
			#print(d)
		elif 'TITLE' in line:
			pass
		elif 'RTINSECONDS' in line:
			pass
		elif 'PEPMASS' in line:
			pass
		elif 'CHARGE' in line:
			pass
		else:
			if 'END IONS' not in line:
				#(firstWord, rest) = line.split(maxsplit=1)
				firstword = line.split(' ')
				key_value = firstword[0]
				d[key] += key_value + ' '
				#print(d)
			else:
				pass
		#print(d)
	
		line =fh.readline()
		#print(d)
		

	store_values = csv.writer(csvFile, delimiter=',')
	#store_values.writerow(["Peptide","Peptide Mass","Protein Id","Protein Description","Protein length","Start pos","end pos","peptide size"])
	store_values.writerow(["SpecFile","SpecID","ScanNum","FragMethod","Precursor","IsotopeError","PrecursorError","Charge","Peptide","Protein","DeNovoScore","MSGFScore","SpecEValue","EValue","QValue","PepQValue","Target_and_decoy","Target_and_decoy_1","base_peptide","Peptide_length","miss_cleavage",'Modifiable_site_cys','Modifiable_site_oxidation','Total_Modifiable_site','base_peptide_without_modification','Scan_id','mass','b_ion_mass','y_ion_mass','Peptide_mass_with_cyscam','Peptide_mass_with_oxidation','Total_mass_with_modification','total_list','b_ion_total_list1'])

	total_mass = ''
	Mass_list = 0
	for m in csv_reader:
		SpecFile=m[0]	
		SpecID=m[1]	
		ScanNum=m[2]	
		FragMethod=m[3]	
		Precursor=m[4]	
		IsotopeError=m[5]	
		PrecursorError=m[6]	
		Charge=m[7]	
		Peptide=m[8]	
		Protein=m[9]	
		#Target_and_Decoy=m[10]	
		#Target_and_Decoy_1=m[11]	
		DeNovoScore=m[10]	
		MSGFScore=m[11]	
		SpecEValue=m[12]
		EValue	=m[13]
		QValue=m[14]	
		PepQValue=m[15]

		#b = m[5]
		#mass = m[4]
		#charge = m[2]
		
		Target_and_decoy=''
		row = 1
		for num in SpecFile:
			#make row of decoy and target
			if 'XXX_' in Protein:
				Target_and_decoy = 'DECOY'
				Target_and_decoy_1 = '0'
			else:
				Target_and_decoy = 'TARGET'
				Target_and_decoy_1 = '1'

			#creat base peptide 
			base_peptide=Peptide[2:-2]
			base_peptide_without_modification1 = ''.join([i for i in base_peptide if not i.isdigit()])
			base_peptide_without_modification2 = base_peptide_without_modification1.replace("+", "")
			base_peptide_without_modification = base_peptide_without_modification2.replace(".", "")
			pep_length= len(base_peptide_without_modification)
			#peptide_list=list(base_peptide)
			#peptide_list1 = peptide_list[::-1]
			#print(peptide_list)
			#calculate Miss Cleavage
			miss_cleavage_K = base_peptide.count("K")
			miss_cleavage_R = base_peptide.count("R")
			miss_cleavage_KP = base_peptide.count(('KP'))
			miss_cleavage_RP = base_peptide.count("RP")
			if base_peptide[-1] == 'K' or  base_peptide[-1] == 'R':
				miss_cleavage = miss_cleavage_K + miss_cleavage_R - miss_cleavage_KP - miss_cleavage_RP - 1
			else:
				miss_cleavage = miss_cleavage_K + miss_cleavage_R - miss_cleavage_KP - miss_cleavage_RP
			
			#calculate modifieble site
			Modifiable_site_cys = base_peptide.count("+57.021")
			Modifiable_site_oxidation = base_peptide.count("+15.995")
			Total_Modifiable_site = Modifiable_site_cys + Modifiable_site_oxidation
			#miss_cleavage = miss_cleavage_K + miss_cleavage_R - miss_cleavage_KP - miss_cleavage_RP - 1
			#miss_cleavage = miss_cleavage_KP
			
			#make colume of only scan id
			Scan_id = SpecID.replace('index=','')

			#mgf file peak list
			str2 = int(Scan_id)
			peptide_mass = d[str2]
			peptide_mass1 = peptide_mass.split()
			#print(peptide_mass1)
			#print(total_list)
			#print(total_list1)

			#Calculate modified site mass and mass of protein
			seq = base_peptide_without_modification
			#sorce of amino acid mass - http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html and here 'c' = cyscam, 'o' = oxidation. 'h' = h2o
			weights = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333, 'c':57.021 , 'o':15.995, 'h' :18.01528}
			weights1 = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333, 'c':57.021 , 'o':15.995, 'h' :18.01528}
			mass = sum(weights[p] for p in seq)
			Mass_list = 0
			total_list = []
			total_list1=[]
			for p in seq:
				mass_list = weights[p]
				Mass_list +=mass_list
				total_list.append(Mass_list)
				#print(total_list)
			
			#b ion calculate
			seq1 = base_peptide.replace('+57.021','c')
			seq1 = seq1.replace('+15.995','o')
			Mass_list1=1.0080
			b_ion_total_list = []
			mass_list1=0
			#mass_list1 = float(mass_list1)
			#Mass_list1 = float(Mass_list1)
			#Mass_list1 = float(Mass_list1)
			for b_ion in seq1:
				if seq1 == 'c':
					mass_list1 = weights1[b_ion]
					Mass_list1 +=mass_list1
					#total_list1.append(Mass_list)
					
				elif seq1 == 'o':
					mass_list1 = weights1[b_ion]
					Mass_list1 +=mass_list1
					#total_list1.append(Mass_list)
					
				else:
					mass_list1 = weights1[b_ion]
					Mass_list1 +=mass_list1
					b_ion_total_list.append(Mass_list1)
					#print(b_ion_total_list)



			#y ion calculate
			seq2 = base_peptide.replace('+57.021','c')
			seq2 = seq1.replace('+15.995','o')
			#seq3 = 'h'
			#seq2 = seq3+seq2
			#print(seq2)
			Mass_list2=b_ion_total_list[-1] 
			Mass_list2 += weights1['h']
			#print(Mass_list2)
			y_ion_total_list = []
			mass_list2=0
			#mass_list1 = float(mass_list1)
			#Mass_list1 = float(Mass_list1)
			#Mass_list1 = float(Mass_list1)
			for y_ion in seq2:
				if seq2 == 'c':
					mass_list2 = weights1[y_ion]
					Mass_list2 -=mass_list2
					#total_list1.append(Mass_list)
					
				elif seq2 == 'o':
					mass_list2 = weights1[y_ion]
					Mass_list2 -=mass_list2
					#total_list1.append(Mass_list)
					
				else:
					mass_list2 = weights1[y_ion]
					y_ion_total_list.append(Mass_list2)
					Mass_list2 =Mass_list2 - mass_list2
					
					#print(Mass_list2)
					
					#print(y_ion_total_list)
					


			#match experimetal peak vs calculated peak

			#common_list = total_list.intersection(peptide_mass)
			#common_list = set(total_list).intersection(*map(set, set(peptide_mass)))
				#common_list = total_list.intersection(d[Scan_id])
			#Total_mass_match = len(common_list)
			

			
			#calculate b_ion and y_ion mass
			b_ion_mass = mass-17
			y_ion_mass = mass+1

			#calculate peptide mass with misscleavage
			Peptide_mass_with_cyscam = Modifiable_site_cys*+57.021+mass
			Peptide_mass_with_oxidation = Modifiable_site_oxidation*+15.995+mass
			Total_mass_with_modification = Modifiable_site_cys*+57.021+Modifiable_site_oxidation*+15.995+mass
			

			
				

			#Scan_id = SpecID.replace('index=','')
			
			store_values.writerow([SpecFile,SpecID,ScanNum,FragMethod,Precursor,IsotopeError,PrecursorError,Charge,Peptide,Protein,DeNovoScore,MSGFScore,SpecEValue,EValue,QValue,PepQValue,Target_and_decoy,Target_and_decoy_1,base_peptide,pep_length,miss_cleavage,Modifiable_site_cys,Modifiable_site_oxidation,Total_Modifiable_site,base_peptide_without_modification,Scan_id,mass,b_ion_mass,y_ion_mass,Peptide_mass_with_cyscam,Peptide_mass_with_oxidation,Total_mass_with_modification,Mass_list,total_list,b_ion_total_list,y_ion_total_list])
			row+=1
			break;