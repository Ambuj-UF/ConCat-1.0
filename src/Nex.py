################################################################################################################
# Handles the rich annotation nexus files                                                                      #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################

import glob
import sys
import csv
from Handler import *
from Functions import *
import timeit



def richNexusCall(runRNA,
                  includeTax,
                  excludeTax,
                  runShannon,
                  calRCVvalue,
                  addTaxName,
                  remTaxName,
                  pipeID,
                  RYcodingCall,
                  spellScan,
                  runBlock,
                  cutOff,
                  rbin,
                  ebin
                  ):
    
    start = timeit.default_timer()
    def transferRNA(file_list):
        file_list = glob.glob("*.nex")
        typeDict = dict()
        RNAfileList = []
        RNAstruc = dict()
        for files in file_list:
            f1 = open(files, 'r')
            flist = f1.readlines()
            parsing = False
            for line in flist:
                if line.rstrip('\n') == "end;":
                    parsing = False
                if parsing:
                    if "Ali_Type = " in line:
                        typeDict[files] = (line.split(" ")[2].rstrip(';\n'))
                    elif "RNA_Type = " in line:
                        if line.split(" ")[2].rstrip(';\n') == 'True':
                            RNAfileList.append(files)
                    elif "RNA_Struc = " in line:
                        RNAstruc[files] = (line.split("=")[1].rstrip(';\n').replace(" ", ""))
                if line.rstrip('\n') == "begin ConCat;":
                    parsing = True

        
        try:
            for files in RNAfileList:
                shutil.copy2(os.path.abspath('') + "/" + files, os.path.abspath('..') + "/" + 'RNAdata')
        except:
            pass
        retVal = [typeDict, RNAstruc]
        return retVal
                

    if pipeID == True:
        usr_inpT = 1
    else:
        usr_inpT = 2

    if usr_inpT == 1:
        
        idDictData = NexusHandler('').managePipes()
        os.chdir("Input/ProcInput")
        file_list = glob.glob("*.nex")
        os.chdir("..")
        if runBlock == True:
            transferRNAret = transferRNA(file_list)
        else:
            transferRNAret = [{},'']
        
        fileTypes = transferRNAret[0]
        os.chdir("ProcInput")
        
        if spellScan == True:
            print "Please wait! Searchng for spelling mistakes \n"
            print "---------------------------------------------------------------------------"
            BaseHandle(2).fuzyName()
            print "---------------------------------------------------------------------------"
        else:
            pass

        try:
            nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
        except:
            inFiles = glob.glob("*.nex")
            for f in inFiles:
                os.remove(f)

            sys.exit("Duplicate alignment files present in Input folder\nProgram Terminated\n")

        combined = Nexus.combine(nexi)
        os.chdir("../..")

    elif usr_inpT == 2:
        os.chdir("Input")
        file_list = glob.glob("*.nex")
        
        if runBlock == True:
            transferRNAret = transferRNA(file_list)
        else:
            transferRNAret = [{},'']
        
        fileTypes = transferRNAret[0]
        
        if spellScan == True:
            print "Please wait! Checking spelling mistakes \n"
            print "---------------------------------------------------------------------------"
            BaseHandle(2).fuzyName()
            print "---------------------------------------------------------------------------"
        else:
            pass

        try:
            nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
        except:
            sys.exit("Duplicate alignment files present in Input folder\nProgram Terminated\n")

        combined = Nexus.combine(nexi)
        os.chdir("..")


    RNAstrucData = transferRNAret[1]

    combinedRet = NexusHandler('filename').NexusHandle(combined,
                                                       usr_inpT,
                                                       RNAstrucData,
                                                       runRNA,
                                                       runShannon,
                                                       includeTax,
                                                       excludeTax
                                                       )
    combined = combinedRet[0]
    remTaxDict = combinedRet[1]
    entropyDict = combinedRet[2]
    for values, data in combined.charsets.items():
        if values.count('.') == 3 and values.split('.')[0] == values.split('.')[2]:
            try:
                del combined.charsets[values]
            except IndexError:
                pass

        elif values.count('.') > 3:
            combined.charsets['.'.join(values.split('.')[1:])] = combined.charsets.pop(values)

    for values, data in combined.taxsets.items():
        if values.count('.') == 3 and values.split('.')[0].split('_')[1] == values.split('.')[2]:
            try:
                del combined.taxsets[values]
            except IndexError:
                pass


    print "Concatenation completed! \n"

    if calRCVvalue == True:
        print "Calculating RCV values for following genes \n"
        newMSA = MultipleSeqAlignment(NexusHandler('fname').combineToRecord(combined))
        rcvDict = dict()

        try:
            for key, val in combined.charsets.items():
                try:
                    if key != 'RNA_Stem' or key != 'RNA_Loop' or key != "'RNA_Stem'" or key != "'RNA_Loop'":
                        if any(fileTypes) == True:
                            if fileTypes[key] == 'Protein':
                                print key
                                try:
                                    rcvDict[key] = ("[ %s ]" % RCVprotCal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                                except ZeroDivisionError:
                                    pass
                                except IndexError:
                                    print "Alignment files not found \n"
                                

                            else:
                                print key
                                try:
                                    rcvDict[key] = ("[ %s ]" % RCVcal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                                except ZeroDivisionError:
                                    pass
                                except IndexError:
                                    print "Alignment files not found \n"

                        else:
                            print key
                            try:
                                rcvDict[key] = ("[ %s ]" % RCVcal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                            except ZeroDivisionError:
                                pass
                            except IndexError:
                                print "Alignment files not found \n"
                except KeyError:
                    continue
        except TypeError:
            print "Type error in RCV Calculation. Skipping this step \n"
            pass

        print "RCV Calculation Completed \n"
        
    else:
        pass
    
    file_name = open('Results1.nex', 'w')
    combined.write_nexus_data(file_name)
    file_name.close()

    def repUPPER(str1, str2):
        return str1.replace(str2, str2.upper())

    with open("Results1.nex", "r") as textfile, open('Results.nex', 'w') as myfile:

        d = textfile.readlines()
        for lines in d:
            if 'charset' in lines or 'taxset' in lines or 'charpartition' in lines:
                if 'taxset Database_IDs' in lines:
                    l1 = lines.split(' = ')[1].split(';')[0]
                    str = ''
                    for i, val in enumerate(l1.split(' ')):
                        dataList = [lval for j, lval in enumerate(val.split('|')) if j>=1]
                        newStr = ''
                        for data in dataList:
                            newStr = newStr + '|' + data

                        dataList = newStr
                        val = val.split('|')[0] + "[" + dataList + "]"
                        str = str + ' ' + val

                    lines = lines.split(' = ')[0] + ' = ' + str + ';\n'
        
                str1 = lines
                str2 = lines.split(' ')[0]
                lines = '\t' + repUPPER(str1, str2)
    
            else:
                pass

            myfile.write(lines)


    if RYcodingCall.isatty() == False:
        RYfiles = dict()
        for lines in RYcodingCall:
            RYfiles[lines.rstrip('\n').split(',')[0]] = (lines.rstrip('\n').split(',')[1])

        fullMSAdata = MultipleSeqAlignment(NexusHandler(2).combineToRecord(combined))
        msaRYlist = []
        for key, val in RYfiles.items():
            msaSendObject = fullMSAdata[:, combined.charsets[key][0]:combined.charsets[key][-1]]
            msaOUT = BaseHandle('string').RYcoding(key, val.lstrip(' '), msaSendObject)
            msaRYlist.append(msaOUT)
        
        msaDataRY = msaRYlist[0]
        
        for i, val in enumerate(msaRYlist):
            if i > 0:
                msaDataRY = msaDataRY + val

        with open('RYoutput.phy', 'w') as fp:
            SeqIO.write(msaDataRY, fp, 'phylip-relaxed')
    
        print "Your output is saved in RYoutput.phy file \n"

    try:
        BaseHandle(2).nexML("Results.nex")
        print "Your NEXML file is stored in Results.xml \n"
    except:
        pass

    def nullTest(val):
        if val is not None:
            return True
        else:
            return False


    if nullTest(addTaxName) == True or nullTest(remTaxName) == True:
        d = dict()
        for row in csv.reader(open('Taxanomy.csv')):
            d['%s' % row[0]] = {'Family': row[1], 'Order': row[2], 'Class': row[3], 'Phylum': row[4], 'Kingdom': row[5]}

    taxDict = dict()

    if nullTest(addTaxName) == True:
        nameList = addTaxName.split(',')
        for lines in nameList:
            tID = lines.rstrip('\n')
            for key in d:
                taxDict[key] = [d[key][inkey] for inkey in d[key] if inkey == tID]
    
            taxDict.pop("Species", [tID])
    
            combined = taxanomyClass(taxDict, combined).addTaxanomy()

        sequences = MultipleSeqAlignment(NexusHandler('fname').combineToRecord(combined))
        fopen = open('ResultsEditedTaxon.nex', 'w')
        SeqIO.write(sequences, fopen, "nexus")
        fopen.close()


    elif nullTest(remTaxName) == True:
        counter = 1
        while counter <= int(remTaxName):
            taxDict = dict()
            combined = taxanomyClass(taxDict, combined).remTaxanomy()

            counter = counter + 1

        sequences = MultipleSeqAlignment(NexusHandler('fname').combineToRecord(combined))
        fopen = open('ResultsEditedTaxon.nex', 'w')
        SeqIO.write(sequences, fopen, "nexus")
        fopen.close()

    if cutOff != None:
        print "Searching fast evolving sites"
        fast_evolv_site = fastEvol(combined, cutOff)
        if fast_evolv_site != []:
            with open('Fast_Evolving_Sites', 'w') as fp:
                for val in fast_evolv_site:
                    fp.write("%s\n" %val[0].split('_')[1])
        else:
            print "No fast Evolving site found"

    os.remove("Results1.nex")

    def two2one(list1):
        listx = []
        for val in list1:
            for inval in val:
                listx.append(inval)
        return listx

    if rbin != None or ebin != None:
        binRetData = binAll(rbin, ebin, combined, rcvDict, entropyDict)

    with open("Combined.nex", 'w') as fp:
        file1 = open("Results.nex", 'r')
        d = file1.readlines()
        list1 = []
        list2 = []
        for lines in d:
            if "CHARSET" in lines:
                list1.append(lines)
            if "CHARSET RNA_" in lines:
                list2.append(lines)
        newList = []
        newlines = [newline for newline in list1 if newline not in list2]
        newList.append(newlines)
        newList.append(list2)

        for i, data in enumerate(newList):
            for j, inData in enumerate(data):
                data[j] = inData.rstrip('\n').lstrip('\t')
            newList[i] = data
    
        newList = two2one(newList)

        if calRCVvalue == True:
            for key, val in rcvDict.items():
                for i, lineVal in enumerate(newList):
                    if key == lineVal.split(' ')[1] or "'" + key + "'" == lineVal.split(' ')[1]:
                        newstrng = "[ RCV Score :" + val + "];"
                        lineVal = " ".join((lineVal.rstrip(';'), newstrng))
                    newList[i] = lineVal

        if runShannon == True:
            for key, val in entropyDict.items():
                for i, lineVal in enumerate(newList):
                    if key == lineVal.split(' ')[1] or "'" + key + "'" == lineVal.split(' ')[1]:
                        newstrng = "[ Entropy : %s ];" %val
                        lineVal = " ".join((lineVal.rstrip(';'), newstrng))
                    newList[i] = lineVal


        for key, val in fileTypes.items():
            for i, lineVal in enumerate(newList):
                if key == lineVal.split(' ')[1] or "'" + key + "'" == lineVal.split(' ')[1]:
                    newstrng = "[ Alignment Type = " + val + "];"
                    lineVal = " ".join((lineVal.rstrip(';'), newstrng))
                newList[i] = lineVal
                    


        listTax = []
        listTax1 = []
        listTax2 = []

        for lines in d:
            if "TAXSET" in lines:
                listTax.append(lines)
            if "TAXSET Database_IDs_" in lines:
                listTax1.append(lines)
            if "TAXSET Missing_" in lines:
                listTax2.append(lines)

        listTAXarrange = []

        newlineTax = [newline for newline in listTax if newline not in listTax1 and newline not in listTax2]

        listTAXarrange.append(newlineTax)
        listTAXarrange.append(listTax1)
        listTAXarrange.append(listTax2)

        for i, data in enumerate(listTAXarrange):
            for j, inData in enumerate(data):
                data[j] = inData.rstrip('\n').lstrip('\t')
            listTAXarrange[i] = data

        listTAXarrange = two2one(listTAXarrange)


        counter = 0
        counterTax = 0
        for lines in d:
            if "CHARSET" in lines:
                fp.write("\t%s\n" % newList[counter])
                counter = counter + 1
            elif "TAXSET" in lines:
                fp.write("\t%s\n" % listTAXarrange[counterTax])
                counterTax = counterTax + 1
            else:
                fp.write(lines)
    
    
        try:
            fp.write("\nBegin ConCat Set; \n\t%s = %s;\n" % ("REMOVED_TAX", remTaxDict['REMOVED_TAX']))
            flag = True
        except:
            flag = False
            pass
        

        def reverseDict(d):
            d2 = dict()
            for key, val in d.items():
                for inkey, inval in val.items():
                    d2[inkey] = []

            for key, val in d.items():
                for inkey, inval in val.items():
                    d2[inkey].append(key + " : " + inval)
        
            return d2

        if usr_inpT == 1:
            idDict = reverseDict(idDictData)

            for key, val in idDict.items():
                if key not in combined.taxlabels:
                    fp.write("\t%s = %s;\n" % ("Database_IDs_" + key, val))

        if cutOff != None and fast_evolv_site != [] and flag == True:
            fp.write("\t%s = %s;\n" % ("Fast_Evolving_Sites", fast_evolv_site))

        elif cutOff != None and fast_evolv_site != [] and flag == False:
            fp.write("\nBegin ConCat Set;\n")
            fp.write("\t%s = %s;\n" % ("Fast_Evolving_Sites", fast_evolv_site))
            fp.write("end;\n")
    
        if flag == True:
            fp.write("end;")

        if rbin != None or ebin != None:
            fp.write("begin ConCat_Bin;\n")
            if rbin != None:
                for val in binRetData[0]:
                    fp.write("\t%s;\n" %val)
            if ebin != None:
                for val in binRetData[1]:
                    fp.write("\t%s;\n" %val)
            fp.write("end;\n")

        fp.close()

    os.remove("Results.nex")

    if any(fileTypes) == True:
        f = open("Partition.txt", 'w')
        for key, value in fileTypes.iteritems():
            if key in combined.charsets:
                f.write('%s, %s = %s-%s\n' % (value, key.split('.')[0], combined.charsets[key][0], combined.charsets[key][-1]))
        f.close()

    # Write IDs to csv file
    if usr_inpT == 1:

        f = open('DatabaseID.csv', 'w')
        writer = csv.writer(f, delimiter = ',')

        writer.writerow(['Species'] + [values.split(':')[0].split('.')[0] for values in idDict[idDict.keys()[1]]])
        for key, value in idDict.iteritems():
            writer.writerow([key] + [values.split(':')[1] for values in value])

        f.close()

    try:
        os.remove("RNAConsensus.txt")
    except OSError:
        pass

    try:
        os.remove("alirna.ps")
    except OSError:
        pass

    try:
        inFiles = glob.glob("*.aln")
        for f in inFiles:
            os.remove(f)
    except:
        pass


    os.chdir("Input/ProcInput")
    inFiles = glob.glob("*.nex")
    for f in inFiles:
        os.remove(f)
    
    print "Your final concatenated result is saved in Combined.nex \n Have a nice day!!"

    os.chdir("../../RNAdata")
    files = glob.glob("*.*")
    try:
        for filename in files:
            if filename.split(".")[1] == 'py' or filename.split(".")[1] == 'pyc' or filename.split(".")[1] == 'md':
                continue
            else:
                try:
                    os.remove(filename)
                except:
                    pass
    except:
        pass
    stop = timeit.default_timer()
    print "Your final concatenated alignment is saved in Combined.nex \n Have a nice day!!"
    print stop - start

















