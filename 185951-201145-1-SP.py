import sys
import math
import os
from scipy import stats

#def ttest(set1, set2):
#    tt = (find_mean(set1) - find_mean(set2)) / np.sqrt(sv/float(n))  # t-statistic for mean
#    pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)

def find_mean(values):
    sum = 0
    for value in values:
        sum += value

    if len(values) != 0:
        return sum/float(len(values))
    else:
        return None

def find_stdev(values):
    mean = find_mean(values)
    sum = 0
    for value in values:
        sum += (value - mean)**2

    if len(values) != 0:
        return math.sqrt(sum/float(len(values)))
    else:
        return None

def get_comparison_dict():
    comparisonDict = {}

    comparisonFile = open("comparisons.txt")
    for line in comparisonFile:
        row = line.split("\t")
        comparisonDict[row[0].strip()] = {"HD" : row[1].strip(),
                                  "control" : row[2].strip()}

    return comparisonDict

def analyze(path, seekGenes, allGenesDict = {}, heatmapGenesDict = {}, batchDict = {}, huntingtonDict = {}):
    fileHandle = open(path)
    
    fileName = path[path.rfind('/') + 1 : path.rfind('_')]
    batch = fileName.replace('GDS', '')

    subsetItems = []
    subsetItemDict = {}
    comparisonDict = get_comparison_dict()

    geneDict = {}

    for gene in seekGenes:
        geneDict[gene] = []

    for line in fileHandle:
        if line.startswith("^DATASET"):
            dataset = line.strip().split(" = ")[1]

        if line.startswith("!dataset_platform_organism"):
            organism = line.strip().split(" = ")[1]

        if "^SUBSET" in line:
            subsetItems.append(subsetItemDict)
            subsetItemDict = {}

        if line.startswith("!subset"):
            split = line.split(" = ")
            subsetItemDict[split[0].strip()] = split[1].strip()

        if line.startswith("!dataset_table_begin"):
            subsetItems.append(subsetItemDict)
            break

    del subsetItems[0]
    subsetDict = {}
    subsetIndexDict = {}

    for line in fileHandle:
        if line.startswith("ID_REF"):
            headerLine = line
            headerRow = headerLine.rstrip().split("\t")
            break

    for item in subsetItems:
        #print "{1}: {0}".format(item["!subset_sample_id"], item["!subset_description"])
        subsetList = item["!subset_sample_id"].strip().split(",")
        subsetDict[item["!subset_description"]] = subsetList

        subsetIndexList = []
        for subject in subsetList:
            subsetIndexList.append(headerRow.index(subject))

        subsetIndexDict[item["!subset_description"]] = subsetIndexList

    while not line.startswith("!dataset_table_end"):
        line = fileHandle.next()
        row = line.rstrip().split("\t")
        if len(row) <= headerRow.index("Gene symbol"):
            continue

        geneName = row[headerRow.index("IDENTIFIER")].strip().upper()
        if geneName.strip() == "":
            continue
        
        means = []
        deviations = []
        counts = []

        geneValueDict = {}

        for subset in subsetIndexDict:
            geneValues = []        
            indices = subsetIndexDict[subset]

            for index in indices:
                value = row[index]
                if value != "null":
                    geneValues.append(float(row[index]))
            
            means.append(find_mean(geneValues))
            deviations.append(find_stdev(geneValues))
            counts.append(len(geneValues))

            geneValueDict[subset] = geneValues

        if dataset in comparisonDict:
            normalSubset = comparisonDict[dataset]["control"]
            huntingtonSubset = comparisonDict[dataset]["HD"]
            
            
			
            if len(geneValueDict[normalSubset]) == 0 or len(geneValueDict[huntingtonSubset]) == 0:
                continue

            tTestResult = stats.ttest_ind(geneValueDict[huntingtonSubset], geneValueDict[normalSubset])


            normalMean = find_mean(geneValueDict[normalSubset])
            huntingtonMean = find_mean(geneValueDict[huntingtonSubset])

            # Include data into heatmapGenesDict
            if not geneName in heatmapGenesDict:
                heatmapGenesDict[geneName] = {}
            for sampleName in subsetDict[normalSubset]:
                batchDict[sampleName] = batch
                huntingtonDict[sampleName] = 0
                if not sampleName in heatmapGenesDict[geneName]:
                    heatmapGenesDict[geneName][sampleName] = []
                value = row[headerRow.index(sampleName)]
                if value and value != "null":
                    heatmapGenesDict[geneName][sampleName].append(float(value))
                else:
                    continue
            for sampleName in subsetDict[huntingtonSubset]:
                batchDict[sampleName] = batch
                huntingtonDict[sampleName] = 1
                if not sampleName in heatmapGenesDict[geneName]:
                    heatmapGenesDict[geneName][sampleName] = []
                value = row[headerRow.index(sampleName)]
                if value and value != "null":
                    heatmapGenesDict[geneName][sampleName].append(float(value))
                else:
                    continue
                    
            if not geneName in allGenesDict:
                allGenesDict[geneName] = {}
               
            #change = (huntingtonMean - normalMean) / normalMean
            change = math.log(huntingtonMean/normalMean)
            pval = tTestResult[1]
            if change < 0:
                pval = -1 * pval
            allGenesDict[geneName][dataset] = [change, pval]

        else:
            tTestResult = ["", ""]
            normalSubset = 0
            huntingtonSubset = 0

        if geneName in seekGenes:
            geneDict[geneName] = {}
            geneDict[geneName]["subsets"] = subsetIndexDict.keys()
            geneDict[geneName]["study"] = dataset
            geneDict[geneName]["organism"] = organism
            geneDict[geneName]["means"] = means
            geneDict[geneName]["deviations"] = deviations
            geneDict[geneName]["counts"] = counts
            geneDict[geneName]["t.statistic"] = tTestResult[0]
            geneDict[geneName]["p.value"] = tTestResult[1]

    if normalSubset:
        huntingtonCount = len(geneValueDict[huntingtonSubset])
        controlCount = len(geneValueDict[normalSubset])
    else:
        huntingtonCount = 0
        controlCount = 0
    return [geneDict, allGenesDict, huntingtonCount, controlCount, heatmapGenesDict, batchDict, huntingtonDict]


def table_to_tsv(table):
    outputString = ""
    for row in table:
        outputString += "\t".join(row) + "\n"

    return outputString

def all_genes_dict_to_table(allGenesDict):
    pval_table = []
    change_table = []

    studyList = allGenesDict[allGenesDict.keys()[0]].keys()
    headerRow = ["gene"]
    for study in studyList:
        headerRow.append(study)

    for gene in allGenesDict:
        pval_row = [gene]
        change_row = [gene]
        for col in headerRow:
            pval_row.append("")
            change_row.append("")

        for study in allGenesDict[gene]:
            if not study in studyList:
                studyList.append(study)
                headerRow.append(study)
                pval_row.append("")
                change_row.append("")

            values = allGenesDict[gene][study]        
            change = values[0]
            pval = values[1]

            colIndex = headerRow.index(study)
            pval_row[colIndex] = str(pval)
            change_row[colIndex] = str(change)

        pval_table.append(pval_row)
        change_table.append(change_row)

    change_table.insert(0, headerRow)
    pval_table.insert(0, headerRow)

    changeTSV = table_to_tsv(change_table)
    pvalTSV = table_to_tsv(pval_table)

    return (changeTSV, pvalTSV)


def heatmap_dict_to_table(heatmapDict, batchDict, huntingtonDict):
    table = []    

    #header row
    headerRow = ["Gene"]
    batchRow = ["Batch"]
    huntingtonRow = ["HD"]

    for gene in heatmapDict:
        row = [gene]
        geneDict = heatmapDict[gene]
        for sample in geneDict:
            if not sample in headerRow:
                headerRow.append(sample)
                batchRow.append(str(batchDict[sample]))
                huntingtonRow.append(str(huntingtonDict[sample]))
            values = geneDict[sample]
            
            if len(values) == 0:
                continue
            
            sum = 0
            for value in values:
                sum += value
            value = str(sum / len(values))
            
            while len(row) < len(headerRow):
                row.append("")
            
            row[headerRow.index(sample)] = value
            
        table.append(row)
    
    table.insert(0, headerRow)
    table.insert(1, batchRow)
    table.insert(2, huntingtonRow)
    #print table
            
    return table_to_tsv(table)
   
def to_table(geneDictList, seekGenes):
    outputTable = []
		
    # Determine how many groups
    groupCount = 0
    for geneDict in geneDictList:
        for gene in geneDict:
            if not geneDict[gene]:
                print "Gene {0} was not found in one of the studies".format(gene)
                continue
            if len(geneDict[gene]["subsets"]) > groupCount:
                groupCount = len(geneDict[gene]["subsets"])

    geneOutputTableDict = {}
    for gene in seekGenes:
        if not gene in geneOutputTableDict:
            geneOutputTableDict[gene] = [["", "", gene]]

        geneOutputHeaderRow = []
        geneOutputHeaderRow.append("Study ID")
        geneOutputHeaderRow.append("Organism")
        geneOutputHeaderRow.append("t-statistic")
        geneOutputHeaderRow.append("p value")
        geneOutputHeaderRow.append("HD")
        geneOutputHeaderRow.append("hdMean")
        geneOutputHeaderRow.append("hdCount")
        geneOutputHeaderRow.append("hdSD")
        geneOutputHeaderRow.append("Control")
        geneOutputHeaderRow.append("ctrMean")
        geneOutputHeaderRow.append("ctrCount")
        geneOutputHeaderRow.append("ctrSD")
        geneOutputHeaderRow.append("change")

        geneOutputTableDict[gene].append(geneOutputHeaderRow)

    comparisonDict = get_comparison_dict()
    
    for gene in geneOutputTableDict:
        for geneDict in geneDictList:
            if not geneDict[gene]:
                continue
            row = [geneDict[gene]["study"], geneDict[gene]["organism"], str(geneDict[gene]["t.statistic"]), str(geneDict[gene]["p.value"])]
            
            for i in range(len(geneDict[gene]["subsets"])):
                subset = geneDict[gene]["subsets"][i]
                if subset == comparisonDict[ geneDict[gene]["study"] ]["control"]:
                    controlName = subset
                    controlMean = str(geneDict[gene]["means"][i])
                    controlCount = str(geneDict[gene]["counts"][i])
                    controlSD = str(geneDict[gene]["deviations"][i])
                if subset == comparisonDict[ geneDict[gene]["study"] ]["HD"]:
                    hdName = subset
                    hdMean = str(geneDict[gene]["means"][i])
                    hdCount = str(geneDict[gene]["counts"][i])
                    hdSD = str(geneDict[gene]["deviations"][i])
                
            row.append(hdName)
            row.append(hdMean)
            row.append(hdCount)
            row.append(hdSD)
            row.append(controlName)
            row.append(controlMean)
            row.append(controlCount)
            row.append(controlSD)
            row.append( str((float(hdMean) - float(controlMean)) / float(controlMean)))

            geneOutputTableDict[gene].append(row)

    for gene in geneOutputTableDict:
        for row in geneOutputTableDict[gene]:
            outputTable.append(row)
        outputTable.append([])

    tsv = table_to_tsv(outputTable)
    return tsv

#argv = sys.argv
#path = argv[1]
path = "C:/Users/Eric/Google Drive/URO/2014/REX Literature/expression/files"

#seekGenes = ["NM_005612", "NM_170735"]
seekGenes = ["REST", "BDNF", 'NTRK1', 'NTRK2', 'NTRK3', 'AKT1', 'AKT2', 'AKT3', 'PI3K', 'FASL', 'NGFR', 'PTEN']

geneDictList = []
allGeneDictHash = {}

if not os.path.isdir(path):
    analysis = analyze(path, seekGenes)
    geneDictList.append(analysis[0])
    allGeneDictList.append(analysis[1])

else:
    allGenesDict = {}
    heatmapGenesDict = {}
    batchDict = {}
    huntingtonDict = {}
    hdCount = 0
    ctrCount = 0
    fileList = []
    for filePath in os.listdir(path):
        if filePath.lower().endswith("soft"):
            fileList.append( path[path.rfind('/') + 1 : path.rfind('_')] )
            
            print "RUNNING: {0}".format(filePath)
            analysis = analyze(path + "/" + filePath, seekGenes, allGenesDict, heatmapGenesDict, batchDict, huntingtonDict)
            geneDictList.append(analysis[0])
            allGenesDict = analysis[1]
            heatmapGenesDict = analysis[4]
            batchDict = analysis[5]
            huntingtonDict = analysis[6]
            print "HD: {0}, Normal: {1}".format(analysis[2], analysis[3])

            hdCount += analysis[2]
            ctrCount += analysis[3]
    print "HD: {0}".format(hdCount)
    print "Control: {0}".format(ctrCount)
    
pValueTable = []
pValHeader = ["Gene", "Study", "p value", "HD Mean", "Control Mean", "Percent Difference"]
pValueTable.append(pValHeader)

comparisonDict = get_comparison_dict()

#print geneDictList

#print heatmapGenesDict
table = to_table(geneDictList, seekGenes)
allGenesTables = all_genes_dict_to_table(allGenesDict)
heatmapTable = heatmap_dict_to_table(heatmapGenesDict, batchDict, huntingtonDict)

outputHandle = open("output.txt", "w")
outputHandle.write(table)
outputHandle.close()

outputHandle = open("output_genechanges.txt", "w")
outputHandle.write(allGenesTables[0])
outputHandle.close()

outputHandle = open("output_genepval.txt", "w")
outputHandle.write(allGenesTables[1])
outputHandle.close()

outputHandle = open("output_heatmap.txt", "w")
outputHandle.write(heatmapTable)
outputHandle.close()
