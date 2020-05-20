import pysam
import re
import sys

regex = re.compile('[0-9]+N[0-9]+D[0-9]+N') 
sj = ''
infile  = pysam.AlignmentFile(sys.argv[1], mode="rb")#, threads=4)
outfile = pysam.AlignmentFile(sys.argv[2], mode="wb", template=infile)#, threads=4) #
for s in infile:
    s2 = s
    if s2.has_tag('MQ'):
        MQtag = s2.get_tag('MQ')
        if MQtag==255:
            s2.set_tag( 'MQ', 60 )
            s2.mapping_quality=60
    if s2.has_tag('MC'):#try:
        ctmp = s2.get_tag('MC')
        rctmp = regex.search(ctmp)
        while rctmp!=None:
            ctmpshort = ctmp[rctmp.start():rctmp.end()]
            indsshort  = [i for i in range(len(ctmpshort)) if ctmpshort[i].isalpha()]
            ndn = int( ctmpshort[0:indsshort[0]] ) + int( ctmpshort[(indsshort[0]+1):indsshort[1]]) + int( ctmpshort[(indsshort[1]+1):indsshort[2]] )
            ctmp = ctmp[0:rctmp.start()]+str(ndn)+"N"+ctmp[rctmp.end():len(ctmp)]
            rctmp = regex.search(ctmp)
        s2.set_tag( 'MC', ctmp )
    outfile.write(s2)
infile.close()
outfile.close()