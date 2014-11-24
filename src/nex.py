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
import time
import platform
from handler import *
from functions import *
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
                  ebin,
                  GC,
                  gcbin,
                  pbin,
                  usrGCbin
                  ):
    
    print("INFO                |    %s|    ------------- ConCat v1.0 ------------------" %(time.strftime("%c")))
    print("INFO                |    %s|    -----Copyright (C) {2014} Ambuj Kumar-------" %(time.strftime("%c")))
    print("INFO                |    %s|    ----------Kimball-Braun Lab group-----------" %(time.strftime("%c")))
    print("INFO                |    %s|    -----You are using python version %s-----" %(time.strftime("%c"), platform.python_version()))
    #print("Start ConCat-build |    %s|    %s|    %s" %(time.strftime("%c"), time.strftime("%H:%M:%S"), sys.version))
    
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
                if "begin ConCat;" in line:
                    parsing = True

        try:
            for files in RNAfileList:
                shutil.copy2(os.path.abspath('') + "/" + files, os.path.abspath('..') + "/" + 'RNAdata')
        except:
            pass
        retVal = [typeDict, RNAstruc]
        return retVal

    print("Start ConCat-build  |    %s|    -------------- Loading Files ---------------" %(time.strftime("%c")))

    if pipeID == True:
        usr_inpT = 1
    else:
        usr_inpT = 2

    if usr_inpT == 1:
        print("Running ConCat-build|    %s|    --------- Extracting database ID's ---------" %(time.strftime("%c")))
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
            print("Running ConCat-build|    %s|    Searchng for spelling mistakes" %(time.strftime("%c")))
            #print("Please wait! Searchng for spelling mistakes \n")
            #print("---------------------------------------------------------------------------")
            BaseHandle(2).fuzyName()
            #print("---------------------------------------------------------------------------")
            print("Running ConCat-build|    %s|    Spell check done" %(time.strftime("%c")))
        else:
            pass

        try:
            nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
        except:
            for filename in file_list:
                dataF = open(filename, 'r').readlines()
                os.remove(filename)
                fp = open(filename, 'w')
                flagF = False
                for lines in dataF:
                    if 'BEGIN MacClade;' in lines:
                        flagF = True
                    if flagF == False:
                        fp.write('%s', lines)
                fp.close()
            
            nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
    
        print("Running ConCat-build|    %s|    ----------Initiating concatenation----------" %(time.strftime("%c")))
        combined = Nexus.combine(nexi)
        print("Running ConCat-build|    %s|    ---------- Conatenation completed ----------" %(time.strftime("%c")))
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
            print("Running ConCat-build|    %s|    Searchng for spelling mistakes" %(time.strftime("%c")))
            #print("---------------------------------------------------------------------------")
            BaseHandle(2).fuzyName()
            #print("---------------------------------------------------------------------------")
            print("Running ConCat-build|    %s|    Spell check done" %(time.strftime("%c")))
        else:
            pass

        try:
            nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
        except:
            for filename in file_list:
                dataF = open(filename, 'r').readlines()
                os.remove(filename)
                with open(filename, 'w') as fp:
                    flagF = False
                    for lines in dataF:
                        if 'BEGIN MacClade;' in lines:
                            flagF = True
                        if flagF == False:
                            fp.write('%s' %lines)
            

            nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]

        print("Running ConCat-build|    %s|    -----Initiating concatenation-----" %(time.strftime("%c")))
        combined = Nexus.combine(nexi)
        print("Running ConCat-build|    %s|    ----- Conatenation completed -----" %(time.strftime("%c")))

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


    def _rcvPrint(counterRCV, totalLength):
        if float(counterRCV)/totalLength*100 < 10:
            print("RCV Calculation     |    %s|    %.2f percent completed             |    %s" %(time.strftime("%c"), float(counterRCV)/totalLength*100, key))
        elif 10 <= float(counterRCV)/totalLength*100 < 100:
            print("RCV Calculation     |    %s|   %.2f percent completed             |    %s" %(time.strftime("%c"), float(counterRCV)/totalLength*100, key))
        else:
            print("RCV Calculation     |    %s|  %.2f percent completed             |    %s" %(time.strftime("%c"), float(counterRCV)/totalLength*100, key))


    if calRCVvalue == True:
        
        print("RCV Calculation     |    %s|    ---------- Initiating RCV Calculation ----------" %(time.strftime("%c")))
        
        newMSA = MultipleSeqAlignment(NexusHandler('fname').combineToRecord(combined))
        rcvDict = dict()
        
        totalLength = len(combined.charsets)
        counterRCV = 0

        try:
            for key, val in combined.charsets.items():
                try:
                    if key != 'RNA_Stem' and key != 'RNA_Loop' and key != "'RNA_Stem'" and key != "'RNA_Loop'":
                        if any(fileTypes) == True:
                            try:
                                if fileTypes[key] == 'Protein':
                                    _rcvPrint(counterRCV, totalLength)
                                    
                                    
                                    try:
                                        rcvDict[key] = ("[ %s ]" % RCVprotCal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                                    except ZeroDivisionError:
                                        pass
                                    except IndexError:
                                        print("WARNING             |    %s|Skipping RCV calculation for above file|    %s alignment not found" %(time.strftime("%c"), key))
                                        #print("Alignment files not found \n")
                                        
                                    counterRCV = counterRCV + 1
                                

                                else:
                                    _rcvPrint(counterRCV, totalLength)
                                    #print("RCV Calculation     |    %s|    %s" %(time.strftime("%c"), key))
                                    
                                    try:
                                        rcvDict[key] = ("[ %s ]" % RCVcal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                                    except ZeroDivisionError:
                                        pass
                                    except IndexError:
                                        print("WARNING             |    %s|Skipping RCV calculation for above file|    %s alignment not found" %(time.strftime("%c"), key))
                                        #print("Alignment files not found \n")
                                        
                                    counterRCV = counterRCV + 1
                    
                            except KeyError:
                                
                                _rcvPrint(counterRCV, totalLength)
                                #print("Setting %s alignment type to default alignment type 'DNA'" %key)
                                
                                print("RCV Calculation     |    %s|    %s" %(time.strftime("%c"), key))
                                try:
                                    rcvDict[key] = ("[ %s ]" % RCVcal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                                except ZeroDivisionError:
                                    pass
                                except IndexError:
                                    print("WARNING             |    %s|Skipping RCV calculation for above file|    %s alignment not found" %(time.strftime("%c"), key))
                                    #print("Alignment files not found \n")
                                    
                                counterRCV = counterRCV + 1


                        else:
                            _rcvPrint(counterRCV, totalLength)
                            
                            try:
                                rcvDict[key] = ("[ %s ]" % RCVcal(newMSA[:, combined.charsets[key][0]:combined.charsets[key][-1]]))
                            except ZeroDivisionError:
                                pass
                            except IndexError:
                                print("WARNING             |    %s|Skipping RCV calculation for above file|    %s alignment not found" %(time.strftime("%c"), key))
                                #print("Alignment files not found \n")
                                
                            counterRCV = counterRCV + 1
            
                except KeyError:
                    continue
        except TypeError:
            print("WARNING             |    %s|    Error in RCV Calculation. Skipping this step" %(time.strftime("%c")))
            #print("Type error in RCV Calculation. Skipping this step \n")
            pass
        
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
                        val = val.split('|')[0] + "[" + dataList.lstrip("'").rstrip("'") + "]"
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
        RYfilesImport = RYcodingCall
        for lines in RYfilesImport:
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
    

        #print("Your output is saved in RYoutput.phy file \n")
        print("RY Coding           |    %s|    Your output is saved in RYoutput.phy file" %(time.strftime("%c")))

    try:
        BaseHandle(2).nexML("Results.nex")
        print("Writing Nexml       |    %s|    Your NEXML file is stored in 'Results.xml' file" %(time.strftime("%c")))
        #print("Your NEXML file is stored in 'Results.xml' file \n")
    except:
        pass

    def nullTest(val):
        if val is not None:
            return True
        else:
            return False


    if nullTest(addTaxName) == True:
        
        d = dict()
        
        print("Editing taxa names  |    %s|    Extracting taxanomy data from Taxanomy.csv file" %(time.strftime("%c")))
        #print("Extracting taxanomy data from Taxanomy.csv file\n")
        try:
            for row in csv.reader(open("Taxanomy.csv", 'rU')):
                d['%s' % row[0]] = {'Family': row[1], 'Order': row[2], 'Class': row[3], 'Phylum': row[4], 'Kingdom': row[5]}
        except:
            raise IOError("Taxanomy.csv file not found in ConCat home directory")

        taxDict = dict()

        nameList = addTaxName.split('-')

        for lines in nameList:
            tID = lines.rstrip('\n')
            for key in d:
                taxDict[key] = [d[key][inkey] for inkey in d[key] if inkey == tID]
    
            if "Species" in list(taxDict.keys()):
                taxDict.pop("Species", [tID])
            elif "species" in list(taxDict.keys()):
                taxDict.pop("species", [tID])
            else:
                raise KeyError("Species header not found in Taxanomy.csv file. Check the first coloumn header.")
    
            combined = taxanomyClass(taxDict, combined).addTaxanomy()

        sequences = MultipleSeqAlignment(NexusHandler('fname').combineToRecord(combined))
        fopen = open("ResultsEditedTaxon.nex", 'w')
        SeqIO.write(sequences, fopen, "nexus")
        fopen.close()


    elif nullTest(remTaxName) == True:
        counter = 1
        while counter <= int(remTaxName):
            taxDict = dict()
            combined = taxanomyClass(taxDict, combined).remTaxanomy()

            counter = counter + 1

        sequences = MultipleSeqAlignment(NexusHandler('fname').combineToRecord(combined))
        fopen = open("ResultsEditedTaxon.nex", 'w')
        SeqIO.write(sequences, fopen, "nexus")
        fopen.close()

    if cutOff != None:
        
        #print("Searching fast evolving sites")
        fast_evolv_site = fastEvol(combined, cutOff)
        if fast_evolv_site != []:
            with open('Fast_Evolving_Sites', 'w') as fp:
                for val in fast_evolv_site:
                    fp.write("%s\n" %val[0].split('_')[1])
        else:
            print("Fast Evolving Sites |    %s|    No fast Evolving site found" %(time.strftime("%c")))
            #print("No fast Evolving site found")

    if GC == True:
        gcDict = GCcontent(combined)

    if usrGCbin != None:
        usrGCdict = gcUserBin(int(usrGCbin), gcDict)

    os.remove("Results1.nex")

    def two2one(list1):
        listx = []
        for val in list1:
            for inval in val:
                listx.append(inval)
        return listx

    if rbin != None or ebin != None or gcbin != None:
        if calRCVvalue == False:
            rcvDict = dict()
        if runShannon == False:
            entropyDict == dict()
        if GC == False:
            gcDict = dict()
        
        binRetData = binAll(rbin, ebin, combined, rcvDict, entropyDict, gcDict, gcbin)

    elif pbin == True:
        if calRCVvalue == False:
            rcvDict = dict()
        if runShannon == False:
            entropyDict == dict()
        if GC == False:
            gcDict = dict()

        binData = binPercent(rcvDict, entropyDict, gcDict, combined, calRCVvalue, runShannon, GC)

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
                        newstrng = "[ RCV Score %s ];" %val
                        lineVal = " ".join((lineVal.rstrip(';'), newstrng))
                    newList[i] = lineVal

        if runShannon == True:
            for key, val in entropyDict.items():
                for i, lineVal in enumerate(newList):
                    if key == lineVal.split(' ')[1] or "'" + key + "'" == lineVal.split(' ')[1]:
                        newstrng = "[ Entropy : %s ];" %val
                        lineVal = " ".join((lineVal.rstrip(';'), newstrng))
                    newList[i] = lineVal

        if GC == True:
            for key, val in gcDict.items():
                for i, lineVal in enumerate(newList):
                    if key == lineVal.split(' ')[1] or "'" + key + "'" == lineVal.split(' ')[1]:
                        newstrng = "[ GC content (in percentage) : %s ];" %val
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

        if rbin != None or ebin != None or gcbin != None:
            fp.write("\nbegin ConCat_Bin;\n")
            if rbin != None:
                fp.write("\n\t[Bin RCV data]\n")
                for val in binRetData[0]:
                    fp.write("\t%s;\n" %val)
            if ebin != None:
                fp.write("\n\t[Bin Entropy data]\n")
                for val in binRetData[1]:
                    fp.write("\t%s;\n" %val)
            if gcbin != None:
                fp.write("\n\t[Bin GC data]\n")
                for val in binRetData[2]:
                    fp.write("\t%s;\n" %val)
            fp.write("end;\n")

        elif pbin == True:
            fp.write("\n\nbegin ConCat_Bin;\n")
            if calRCVvalue == True:
                fp.write("\n\t[RCV Bin]\n")
                for key, val in binData[0].items():
                    if key == '0_to_25':
                        fp.write("\n\t[0 to 25th percentile RCV data]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                    elif key == '25_to_75':
                        fp.write("\n\t[25th to 75th percentile RCV data]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                    else:
                        fp.write("\n\t[75th to 100 percentile RCV data]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)


            if runShannon == True:
                fp.write("\n\t[Entropy Bin]\n")
                for key, val in binData[1].items():
                    if key == '0_to_25':
                        fp.write("\n\t[0 to 25th percentile Entropy data]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                    elif key == '25_to_75':
                        fp.write("\n\t[25th to 75th percentile Entropy data]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                    else:
                        fp.write("\n\t[75th to 100 percentile Entropy data]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                                            

            if GC == True:
                fp.write("\n\t[GC Bin]\n")
                for key, val in binData[2].items():
                    if key == '0_to_25':
                        fp.write("\n\t[0 to 25th percentile GC content]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                    elif key == '25_to_75':
                        fp.write("\n\t[25th to 75th percentile GC Content]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)
                    else:
                        fp.write("\n\t[75th to 100 percentile GC content]\n")
                        for inval in val:
                            fp.write("\t%s;\n" %inval)


            if usrGCbin != None:
                fp.write("\n\t[User percentile GC Bin]\n")
                for key, val in usrGCdict.items():
                    fp.write("\t%s = %s\n" %(key, val))

            fp.write("end;\n")

        fp.close()

    os.remove("Results.nex")

    def _checkIndex(listObject, element):
        for i, val in enumerate(listObject):
            if val == element:
                return i

    if any(fileTypes) == True:
        
        print("Creating Partition  |    %s|    Writing alignment partition data to 'Partition.txt' file" %(time.strftime("%c")))
        #print("Writing alignment partition data to 'Partition.txt' file\n")
        f = open("Partition.txt", 'w')
        for key, value in fileTypes.items():
            if key in combined.charsets:
                f.write('%s, %s = %s-%s\n' % (value, key.split('.')[0], combined.charsets[key][0], combined.charsets[key][-1]))
        f.close()

    # Write IDs to csv file
    if usr_inpT == 1:
        
        print("Database ID         |    %s|    Writing database ID's to 'DatabaseID.csv' file" %(time.strftime("%c")))
        #print("Writing database ID's to 'DatabaseID.csv' file\n")
        
        f = open('DatabaseID.csv', 'w')
        writer = csv.writer(f, delimiter = ',')
        
        geneNames = list()
        for key, val in idDict.items():
            for inval in val:
                geneNames.append(inval.split(" :")[0])
    
        geneNames = list(set(geneNames))

        writer.writerow(['Species'] + [x.rstrip(".nex") for x in geneNames])
        for key, value in idDict.items():
            indexData = list()
            for i in range(0, len(geneNames)):
                flag = False
                for inval in value:
                    if inval.split(" : ")[0] == geneNames[i]:
                        indexData.append(inval.split(" : ")[1])
                        flag = True

                if flag == False:
                    indexData.append("NA")
            
            writer.writerow([key] + indexData)

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

    print("Running ConCat-build|    %s|    ------ Processing completed ------" %(time.strftime("%c")))
    print("\nYour final concatenated alignment is saved in Combined.nex \nHave a nice day!!\n")

    if nullTest(addTaxName) == True or nullTest(remTaxName) == True:
        print("ResultsEditedTaxon.nex contains concatenated alignment with edited taxon names\n")

    print("Time elapsed: %s seconds" %(stop - start))

















