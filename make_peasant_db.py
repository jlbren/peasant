## make_peasant_db.py
## Code to make a blast database for peasant annotation
## requires NCBI format ptt files for gene information.
## assumes files are formatted correctly.


import csv
import os
import sys
from Bio import SeqIO
import argparse
import time
import datetime

def msg(name=None):
    return '''make_peasant_db.py -f <filename(s)> -p <filename(s)> -F <filename> -o <database name>'''
parser=argparse.ArgumentParser(usage=msg())
parser.add_argument('-f','--fasta',nargs='+',metavar='<filename>',help='FASTA format nucleotide sequence(s). (space separator)')
parser.add_argument('-p','--ptt',nargs='+',metavar='<filename>',help='NCBI PTT format annotation file entered in the same order as FASTA. (space separator)')
parser.add_argument('-F','--csv_file',action="store",help='CSV file of FASTA format nucleotide sequences, NCBI PTT format file.')
parser.add_argument('-db','--database_name',action="store",help='Name of database')
parser.add_argument('--version',action='version',version='%(prog)s 1.0')
args=parser.parse_args()
fastas=args.fasta
ptts=args.ptt
csvFile=args.csv_file
output_file=args.database_name
if csvFile is not None:
    if fastas is not None and ptts is not None:
        parser.error('Either provide a csv file or values for parameters -f and -p.')
    if not os.path.isfile(csvFile):
        parser.error('ERROR: '+csvFile+' does not exist.')
else:
    if fastas is None or ptts is None:
        parser.error('Both FASTA files and PTT files must be supplied.')
    if len(fastas) != len(ptts):
        parser.error('A PTT file must be supplied for each FASTA file.')
if output_file is None:
    parser.error('A name for the database must be provided.')

log_text='python'
for i in sys.argv:
    log_text+=' '+i

if csvFile is not None:
    fastas=[]
    ptts=[]
    with open(csvFile,'rb') as f:
        reader=csv.reader(f)
        list_of_files=list(reader)
    for i in list_of_files:
        fastas.append(i[0])
        ptts.append(i[1])

new_fasta=open(output_file+'.fasta','w')
for i in range(0,len(fastas)):
    ptt=open(ptts[i],'r')
    ptt.next()  #throw away first line
    ptt.next()  #throw away second line
    labels=ptt.next()
    headers=labels.rstrip('\n').split("\t")

    #get each of the rows
    #each row is saved as a dictionary in a list g
    rows=csv.DictReader(ptt,headers,delimiter='\t')
    g=[]
    for row in rows:
        g.append(row)
    ptt.close()        

    sequences=list(SeqIO.parse(fastas[i],'fasta'))

    #relabel
    #loop through each of the sequences in the fasta file
    for r in sequences:
        new_label=r.name
        position=r.name[r.name.find(':')+1:]
        newposition=position.replace('-','..')  #parse out the position
        if newposition[0]=='c':
            start=newposition[newposition.find('.')+2:]
            end=newposition[1:newposition.find('.')]
            newposition=start+".."+end

        #find the record from the ptt corresponding to this fasta
        found_gene=False
        index_gene=0
        counter=0
        for loop_g in g:
            if loop_g.get('Location')==newposition:
                #found it
                found_gene=True
                index_gene=counter
            counter+=1
        if found_gene==True:
            #format gene_name|synonym|COG|Product_Info|header
            r.name=g[index_gene].get('Gene')+'|'+g[index_gene].get('Synonym')+\
            '|'+g[index_gene].get('COG')+'|'+g[index_gene].get('Product')+'|'+\
            r.name
            r.name=r.name.replace(',',';')
            r.name=r.name.replace(' ','_')
            r.name=r.name.replace('|:','|')
        else:
            r.name='UNKN|UNKN|UNKN|UNKN|'+r.name
            r.name=r.name.replace(',',';')
            r.name=r.name.replace(' ','_')
            r.name=r.name.replace('|:','|')

    for r in sequences:
        new_fasta.write('>'+r.name+'\n')
        strseq=str(r.seq)
        new_fasta.write(strseq+'\n')
new_fasta.close()

blastdbcommand='makeblastdb -in '+output_file+'.fasta -title '+output_file+' -out '+output_file+' -dbtype nucl'
os.system(blastdbcommand)
os.remove(output_file+'.fasta')
ts=time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
log_file=open('make_db_log.txt','w')
log_file.write(st+'\n')
log_file.write('Command: '+log_text)
log_file.close()
