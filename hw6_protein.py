"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    openFile=open(filename,"r")
    text = openFile.read()
    lines=text.splitlines()
    dna=""
    for i in range(len(lines)):
        dna+=lines[i]
    return dna

   

'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    codons=[]
    rna=[]
    dna=dna.replace('T','U')
    for i in range(startIndex,len(dna),3):
        codons.append(dna[i:i+3])
    RNAstartIndex=codons.index("AUG")
    for i in range(RNAstartIndex,len(codons)):
        if codons[i] in ["UAA","UAG","UGA"]:
            rna.append(codons[i])
            break
        rna.append(codons[i])
    return rna


   


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    file=open(filename,"r")
    text=json.load(file)
    codonToAminoAcidDict={}
    for aminoAcid,codonList in text.items():
        for codon in codonList:
            codon=codon.replace('T','U')
            codonToAminoAcidDict[codon]=aminoAcid
    return codonToAminoAcidDict

   


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein=["Start"]
    for i in range(1,len(codons)):
        if codons[i] in ["UAA","UAG","UGA"]:
            protein.append("Stop")
            break
        protein.append(codonD[codons[i]])
    return protein



'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna=readFile(dnaFilename)
    codonDict=makeCodonDictionary(codonFilename)
    proteins=[]
    startValue=0
    while(startValue<len(dna)):
        dnaStr=dna[startValue:]
        dnaStartIndex=dnaStr.find("ATG")
        if (dnaStartIndex<0):
            break
        rnaStrand=dnaToRna(dnaStr, dnaStartIndex)
        protein=generateProtein(rnaStrand,codonDict)     
        proteins.append(protein)
        startValue+=3*len(protein)+dnaStartIndex
    print(len(proteins))
    return proteins

 


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    commonProteinsList=[]
    for list in proteinList1:
        if list in proteinList2:
            commonProteinsList.append(list)
    return commonProteinsList


   


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    aminoAcids=[]
    for protein in proteinList:
        for aminoAcid in protein:
            aminoAcids.append(aminoAcid)
    return aminoAcids


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoAcidDict={}
    for aminoAcid in aaList:
            if aminoAcid in aminoAcidDict:
                aminoAcidDict[aminoAcid]+=1
            else:
                aminoAcidDict[aminoAcid]=1

    return aminoAcidDict

   


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    result=[]
    proteins1=combineProteins(proteinList1)
    proteins2=combineProteins(proteinList2)
    aminoDict1=aminoAcidDictionary(proteins1)
    aminoDict2=aminoAcidDictionary(proteins2)
    count1=len(proteins1)
    count2=len(proteins2)
    for amino in aminoDict1:
        if amino not in aminoDict2:
            aminoDict2[amino]=0.0
        aminoDict1[amino]/=count1
            
    for amino in aminoDict2:
        if amino not in aminoDict1:
            aminoDict1[amino]=0.0
        aminoDict2[amino]/=count2  
            
    
    for aminoAcid in aminoDict1:
        if aminoAcid not in["Start","Stop"]:
                freq1= aminoDict1[aminoAcid]
                freq2= aminoDict2[aminoAcid]
                if abs(freq1-freq2)> cutoff:
                    temp=[]
                    temp.append(aminoAcid)
                    temp.append(freq1)
                    temp.append(freq2)
                    result.append(temp)
    return result

    


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences")
    unique=[]
    for protein in commonalities:
        if protein not in unique:
            unique.append(protein)
        else:continue
    for protein in unique:
        protein.remove("Start")
        protein.remove("Stop")
        if len(protein)==1:
                print(protein[0])
        if len(protein)>1:

            proStr=""
            for i in range(len(protein)-1):
                proStr=proStr+protein[i]+"-"
            proStr=proStr+protein[len(protein)-1]
            print(proStr)             
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for diff in differences:
        print(diff[0] +":"+ str(round(diff[1],2))+"% in Seq1, "+str(round(diff[2],2))+"% in Seq2")       

    return



def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    aminoAcids1=combineProteins(proteinList1)
    aminoAcids2=combineProteins(proteinList2)
    uniqueAminoAcids=[]
    for amino in aminoAcids1+aminoAcids2:
        if amino not in uniqueAminoAcids:
            uniqueAminoAcids.append(amino)
        else: continue
    uniqueAminoAcids.sort()
    return uniqueAminoAcids

   

'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
  aminoAcids=combineProteins(proteinList)
  aminoAcidDict=aminoAcidDictionary(aminoAcids)
  aminoAcidFreq=[]
  for amino in labels:
        if amino in  aminoAcidDict:
            aminoAcidFreq.append(aminoAcidDict[amino]/len(aminoAcids))
        else:
            aminoAcidFreq.append(0)
  print(aminoAcidDict)
  print(aminoAcidFreq)
  return aminoAcidFreq

 


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w = 0.35  # the width of the bars
    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1,edgecolor=edgeList)
    plt.bar(xLabels, freqList2, width=w, align='edge', label=label2,edgecolor=edgeList) 
    title="side by side bar plot"
    plt.xticks(rotation="horizontal")
    plt.legend()
    plt.title(title)
    plt.show()
    return

    


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    colors=[]
    bigDiffList=[]
    for i in range(len(biggestDiffs)):
        bigDiffList.append(biggestDiffs[i][0])
    for label in labels:
        if label in bigDiffList:
            colors.append("black")
        else: colors.append("white")
    print(colors)
    return colors




'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    humanProteins=synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins=synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")
    commonalities=commonProteins(humanProteins,elephantProteins)
    differences=findAminoAcidDifferences(humanProteins,elephantProteins,0.005)
    displayTextResults(commonalities,differences) 
    distinctAminoAcids=makeAminoAcidLabels(humanProteins,elephantProteins)
    humanFreqs=setupChartData(distinctAminoAcids,humanProteins)
    elephantFreqs=setupChartData(distinctAminoAcids,elephantProteins)
    edges=makeEdgeList(distinctAminoAcids,differences)
    createChart(distinctAminoAcids, humanFreqs, "Human", elephantFreqs, "Elephant", edgeList=edges) 
    return





### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
  
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
