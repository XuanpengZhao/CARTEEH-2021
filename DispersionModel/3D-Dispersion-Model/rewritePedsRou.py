
import numpy as np
import queue
# this is the main entry point of this script

f = open("peds.rou.xml", "r")
fwrite = open("peds2.rou.xml", "w")
line = f.readline()
while line:
    if "<walk edges=" in line:
        findQ = False
        source = ""
        for i in line:
            if findQ:
                if i == " " or i == "\"":
                    break
                source += i
            if i == "\"" and not findQ:
                findQ = True
        #print(source)
        line = f.readline()
        if "<ride from=" in line:
            findT = False
            findQ2 = False
            dest = ""
            for i in line:
                if i == "t":
                    findT = True
                if findQ2:
                    if i == "\"":
                        break
                    dest += i
                if findT and i == "\"" and not findQ2:
                    findQ2 = True
            line2write = "        <walk edges=\"" + source + "\" arrivalPos=\"random\"/>\n"
            fwrite.writelines(line2write)
            line2write = "        <ride from=\"" + source + "\" to=\"" + dest + "\" lines=\"taxi\"/>\n"
            fwrite.writelines(line2write)
            #print(dest)
        else:
            print("Didn't find accoding ride")
         
            
    else:
        fwrite.writelines(line)

    line = f.readline()
    
