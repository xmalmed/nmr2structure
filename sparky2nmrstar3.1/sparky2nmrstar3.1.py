#!/usr/bin/python3
# python3
# http://www.bmrb.wisc.edu/software/tablegen/Atom_chem_shift.php?restype=polypeptide&tagcat=Atom_chem_shift
# generate "template.txt" NMR-Star file for primary sequence
# input exported CS list from Sparky
# expecting, first aa in the template is 1

sparkyFile = "sparky.shifts"
templateFile  = "template.txt"
resultFile = "NMRStar-CS.cif"

protein = []

import operator
import re

#------------------------------------------------------------

class AminoAcid:
    
    def __init__(self, number, name):
        self.number = number
        self.name = name
        self.atoms = dict()

#------------------------------------------------------------        

def createModel(templateFile):
    '''
    protein = list of aa instancis
    '''
    protein = []
    head = []
    tail = []
    
    for line in templateFile:
        splitLine = line.split()
        if (len(line) == 0)  or  (line[0] == '#'):
            if len(protein) == 0:
                head.append(line)
            else:
                tail.append(line)
            continue         
        try:
            splitLine[0] = int(splitLine[0])
            number = int(splitLine[2])
            name = splitLine[3]
        except (ValueError, IndexError):
            if len(protein) == 0:
                head.append(line)
            else:
                tail.append(line)
            continue
        
        if len(protein) == 0:
            protein.append(AminoAcid(number,name))
        
        if protein[-1].number != number:
            protein.append(AminoAcid(number,name))
        
        protein[-1].atoms[splitLine[4]] = splitLine
        protein[-1].atoms[splitLine[4]][1] = '1'
        
    return protein, head, tail

#------------------------------------------------------------
    
def printProtein(protein, head, tail, resultFile):
    
    count = 1    
    for line in head:
        resultFile.write(line)
    
    for aa in protein:
        sort = sorted(aa.atoms.items(), key = operator.itemgetter(1))
        for atom in sort:
            if atom[1][7] != '@':
                resultFile.write('{:<5}{:3}{:4}{:4}{:5}{:2}{:3}{:8}{:8}{:3}{:3}{:3}{:3}\n'.format(count+6,atom[1][1],atom[1][2],atom[1][3],atom[1][4],atom[1][5],atom[1][6],atom[1][7],atom[1][8],atom[1][9],atom[1][10],atom[1][11],atom[1][12]))
                count += 1
                
    for line in tail:
        resultFile.write(line)       

#------------------------------------------------------------

def fillIn(protein, sparkyFile):

    rest = []
    for line in sparkyFile:        

        subLine = re.sub(r'^( *[A-Z])', r'\1  ', line)
        splitLine = subLine.split()        
        try:
            res = int(splitLine[1]) - 1
        except (ValueError, IndexError):
            continue
        
        
        
        # out of range of template is ignoring
        try:
            protein[res].number
        except IndexError:
            continue
        
        # expecting, first aa is 1, so position in protein [0]
        if protein[res].number != res + 1:
            raise ValueError('aa is not matching!')

        #specificke upravy...    
        if (splitLine[0] == 'T') and (splitLine[2] == 'MG2'):
            splitLine[2] = 'MG'
        
        if (splitLine[0] == 'I'):
            if splitLine[2] == 'MG2':
                splitLine[2] = 'MG'
            elif splitLine[2] == 'MD1':
                splitLine[2] = 'MD'        
                
        try:    
            if (protein[res].atoms[splitLine[2]][11] == ".") or ( int(protein[res].atoms[splitLine[2]][11]) < int(splitLine[6]) ):
                protein[res].atoms[splitLine[2]][11] = splitLine[6]
                protein[res].atoms[splitLine[2]][7:9] = splitLine[4:6]
                if protein[res].atoms[splitLine[2]][10] == '.':
                    protein[res].atoms[splitLine[2]][10] = '2'
        except KeyError:
            rest.append(splitLine) 
                    
    return protein, rest

#------------------------------------------------------------

def solveRest(protein, rest):
    
    rest2 = []
    for line in rest:
        res = int(line[1]) - 1
        
        if re.match(r'^Q([ABG]1?)$', line[2]):
            k, l, H = '2', '3', 'H'
            vzor = re.match(r'^Q([ABG]1?)$', line[2]).group(1)
            
        elif re.match(r'^Q([DE]2)$', line[2]):
            k, l, H = '1', '2', 'H'
            vzor = re.match(r'^Q([DE]2)$', line[2]).group(1)
        elif re.match(r'^Q(MD)$', line[2]):
            k, l, H = '1', '2', ''
            vzor = re.match(r'^Q(MD)$', line[2]).group(1)   
        elif (re.match(r'^Q([DE])$', line[2])) and (line[0] in 'KRP'):
            k, l, H = '2', '3', 'H'
            vzor = re.match(r'^Q([DE])$', line[2]).group(1)    
        elif (re.match(r'^Q([DE])$', line[2])) and (line[0] in 'FY'): # W
            k, l, H = '1', '2', 'H'
            vzor = re.match(r'^Q([DE])$', line[2]).group(1)  
        elif (re.match(r'^CQ([DE])$', line[2])) and (line[0] in 'FYL'): # CQG Val
            k, l, H = '1', '2', 'C'
            vzor = re.match(r'^CQ([DE])$', line[2]).group(1)  
        else:
            rest2.append(line)
            continue
        
        if (protein[res].atoms[H + vzor + k][7] == '@') and (protein[res].atoms[H + vzor + l][7] == '@'):
            protein[res].atoms[H + vzor + k][7:9] = line[4:6]
            protein[res].atoms[H + vzor + l][7:9] = line[4:6]
            protein[res].atoms[H + vzor + k][10] = '1'
            protein[res].atoms[H + vzor + l][10] = '1'
            protein[res].atoms[H + vzor + k][11] = line[6]
            protein[res].atoms[H + vzor + l][11] = line[6]
        elif (protein[res].atoms[H + vzor + k][7] != '@') and (protein[res].atoms[H + vzor + l][7] != '@'):
            pass
        else:
            rest2.append(line)
                
    return protein, rest2
        

#------------------------------------------------------------
#------------------------------------------------------------

with open(sparkyFile, mode='r', encoding='utf-8') as sparky, open(templateFile, mode='r', encoding='utf-8') as template, open(resultFile, mode='w', encoding='utf-8') as result:
    
    protein, head, tail  = createModel(template)

    protein, rest = fillIn(protein, sparky)
    
    protein, rest = solveRest(protein, rest)
    
    for i in rest:
        print(i)

    printProtein(protein, head, tail, result)    
    
    

