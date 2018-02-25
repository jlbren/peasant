#peasant
import os
import sys
import shutil
import csv
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio import SeqIO
from Bio.SeqUtils import GC
import collections
import time
import datetime
import argparse
import multiprocessing

############ USER NEEDS TO SET THIS PATH ############

database_path='/home/lsb456/peasant/databases'

#####################################################



##################### FUNCTIONALITY TO GET PARAMETERS FROM COMMAND LINE #####################
def msg(name=None):
    return '''peasant.py [assembly options] [filter options] [database] [homology options] -o output_path'''
parser=argparse.ArgumentParser(usage=msg())
groupA=parser.add_mutually_exclusive_group()
groupA.add_argument('-a','--assembler',choices=['spades'],action="store",help='Assembly method.')
groupA.add_argument('-A','--assembled_contigs',action="store",metavar='<filename>',help='Contigs file. Start analysis from an existing assembly.')
groupR=parser.add_mutually_exclusive_group()
groupR.add_argument('-s','--single_reads',action="store",metavar='<filename>',help='Single reads.')
groupR.add_argument('-p','--paired_end_reads', nargs=2, metavar=('<R1 filename>','<R2 filename>'), help='Paired-end reads. List both read files.')
parser.add_argument('-o','--output_path',action="store",metavar='<directory>',help='Directory to store all the resulting files (required)')
parser.add_argument('-t','--num_threads',action="store",metavar='<int>',type=int,help='Number of processors to use.')
groupF = parser.add_argument_group('filter options')
groupF.add_argument('-m','--min_contig_size',action="store",metavar='<int>',type=int,help='Minimum contig size.')
groupF.add_argument('-M','--max_contig_size',action="store",metavar='<int>',type=int,help='Maximum contig size.')
groupF.add_argument('-c','--min_coverage',action="store",metavar='<int>',type=int,help='Minimum coverage.')
groupF.add_argument('-cov','--min_SPAdes_cov',action="store",metavar='<float>',type=float,help='Minimum SPAdes cov value.')
groupDB = parser.add_argument_group('database options')
groupDB.add_argument('-g','--genus',nargs='+',metavar='<filename>',help='BLAST database prefix to use for annotation. To consider more than one database, list genus name (space separator) (required)')
groupHL=parser.add_argument_group('homology thresholds (optional)')
groupHL.add_argument('-q','--qcov',action="store",metavar='<float>',type=float,help='Minimum query coverage to call homologous genes. Default=70.0')
groupHL.add_argument('-i','--pident',action="store",metavar='<float>',type=float,help='Minimum percent identity to call homologous genes. Default=70.0')
groupHL.add_argument('-b','--bitscore',action="store",metavar='<float>',type=float,help='Minimum bitscore to call homologous genes. Default=50.0')
parser.add_argument('--version',action='version',version='%(prog)s 1.0')
args=parser.parse_args()
if args.assembler is None and args.assembled_contigs is None:
    parser.error('Assembler method or assembly file must be specified.')
if args.assembler is not None:
    if args.paired_end_reads is None and args.single_reads is None:
        parser.error('Reads must be provided for assembly.')
if args.output_path is None:
    parser.error('Output path must be provided.')
if args.genus is None:
    parser.error('Database must be provided.')
##################### FUNCTIONALITY TO GET PARAMETERS FROM COMMAND LINE #####################


#PROGRAM FUNCTIONALITY
#CHECK TO MAKE CERTAIN ALL OF THE PATHS ARE SET
def path_check(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

#CHECK ALL PATHS
def check_installs():
    #check sickle
    if path_check('sickle') == None:
        print('Program SICKLE not in path')
        return False

    #check SPAdes
    if path_check('spades.py') == None:
        print('Program SPAdes not in path')
        return False

    #check GLIMMER
    if path_check('glimmer3') == None:
        print('Program GLIMMER not in path')
        return False

    #check tRNAscan-SE
    if path_check('tRNAscan-SE') == None:
        print('Program tRNAscan-SE not in path')
        return False

    #check for blasts (blastx as example)
    if path_check('blastx') == None:
        print('BLAST+ tools are not added to the path')
	return False

    #check BBMAP
    if path_check('bbwrap.sh') == None:
        print('bbmap tools are not added to the path')
	return False


#CHECK IF A FILE EXISTS
def check_for_file(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        return True

#CHECK IF A GROUP OF FILES EXISTS
def check_files(*file_names):
    error_message=''
    for i in file_names:
        if check_for_file(i)==False:
            error_message+='Error: File '+i+' does not exist.\n'
    if len(error_message)==0:
        return True,error_message
    else:
        return False,error_message

#write to log
def log(message,file_name):
    logf=open(file_name,'a')
    logf.write(message+'\n')
    logf.close()

#GET NUMBER OF SEQUENCES IN A FASTA FORMAT FILE
def get_number_of_sequences(contig_file):
    temp_sequences=list(SeqIO.parse(contig_file,'fasta'))
    num_seqs=len(temp_sequences)
    del temp_sequences[:]
    return num_seqs

#CALCULATE N50
def N50(file_name):
    sequences=list(SeqIO.parse(file_name,'fasta'))
    list_of_lens=[]
    sum_of_all=0
    for i in sequences:
        list_of_lens.append(len(i.seq))
        sum_of_all+=len(i.seq)
    list_of_lens.sort(reverse=True)
    halfsies=float(sum_of_all)/2.0
    val=0
    for i in list_of_lens:
        val+=i
        if val>=halfsies:
            return i

#CALCULATE SIZE OF ASSEMBLY
def assembly_size(file_name):
    sequences=list(SeqIO.parse(file_name,'fasta'))
    sum_of_all=0
    for i in sequences:
        sum_of_all+=len(i.seq)
    return sum_of_all

#CALCULATE GC CONTENT
def calc_GC_content(file_name):
    sequences=list(SeqIO.parse(file_name,'fasta'))
    sequence=''
    for i in sequences:
        sequence+=i.seq
    return GC(sequence)

#REFORMATTING OF THE CONTIG HEADERS
def qc_contig_file(contig_file):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    outputfile=open(contig_file,'w')
    for i in sequences:
        description=i.description
        description=re.sub(r'\,','',description)
        description=re.sub(r'\ ','_',description)
        outputfile.write('>'+description+'\n')
        outputfile.write(str(i.seq)+'\n')
    outputfile.close()
    return contig_file

#GENERATE A LIST OF THE CONTIG HEADERS
def get_list_of_contigs(contig_file):
    contigs=[]
    filename=open(contig_file,'r')
    for l in filename:
        if l[0]=='>':
            contigs.append(l[1:].strip())
    filename.close()
    return contigs

#RUN GLIMMER
#make sure the awks are pointing to /usr/bin/awk
def run_glimmer(contig_file,output_path):
    #predict the genes
    glimmer_command='./g3-iterated.csh '+contig_file+' '+output_path+'/temp/orfs'
    print(glimmer_command)
    os.system(glimmer_command)

    #determine which orfs are in which contigs
    contigs=get_list_of_contigs(contig_file)
    contig_label=contigs[0]
    orfs=collections.OrderedDict()
    orfies=open(output_path+'/temp/orfs.coords','r')
    temp_orfs=[i.strip().split() for i in orfies]
    for i in temp_orfs:
        if i[0][0]=='>': #new contig
            contig_label=i[0][1:]
        else:
            orfs[(i[0],contig_label)]=[i[1],i[2],i[3]]
    orfies.close()

    #multi-extract function -- includes the last 3 nucleotides (the stop codon)
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    orf_sequences=[]
    for i,j in orfs: #i is orfID, j is the contig
        for k in sequences:
            if j==k.description:
                start,stop,frame=orfs[(i,j)]
                start=int(start)
                stop=int(stop)
                frame=int(frame)
                if start>stop:
                    if frame<0: #it's reverse complement
                        subseq=k.seq[stop-1:start].reverse_complement()
                    else: #it goes over the origin
                        subseq=k.seq[start-1:]+k.seq[:stop]
                else:
                    if frame<0: #it goes over the origin and reverse complement
                        subseq=(k.seq[stop-1:]+k.seq[:start]).reverse_complement()
                    else:
                        subseq=k.seq[start-1:stop]
                orf_label='>'+j+'_'+i+'_'+str(start)+'_'+str(stop)
                orf_sequences.append([orf_label,subseq])
                break
    orf_output=open(output_path+'/temp/predicted_orfs.fasta','w')
    for i in orf_sequences:
        orf_output.write(str(i[0])+'\n')
        orf_output.write(str(i[1])+'\n')
    orf_output.close()
    return output_path+'/temp/predicted_orfs.fasta'

#RUN READ QC (SICKLE) AND ASSEMBLY (SPADES)
def qc_and_assembly(output_path,assembler,type_reads,*reads):
    if assembler=='spades':
        if type_reads=='pe':
            read1=reads[0][0]
            read2=reads[0][1]         
            if run_sickle(output_path,'pe',read1,read2)==True:
                output_file_label1='trimmed_'+read1[read1.rfind("/")+1:]
                output_file_label2='trimmed_'+read2[read2.rfind("/")+1:]
                if run_spades(output_path,'pe',output_path+'/temp/'+output_file_label1,output_path+'/temp/'+output_file_label2)==True:
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.')
                return False
        else:
            read1=reads[0]
            if run_sickle(output_path,'se',read1)==True:
                output_file_label='trimmed_'+read1[read1.rfind("/")+1:]
                if run_spades(output_path,'se',output_path+'/temp/'+output_file_label)==True:
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.')
                return False


#RUN SICKLE
def run_sickle(output_path,r_type='pe',*read_files):
    sickle_command='sickle '
    if r_type=='pe':    #paired-end reads
        if len(read_files)!=2:
            print('Paired-end reads requires 2 read files entered.')
            return False
        else:
            output_file_label1='trimmed_'+read_files[0][read_files[0].rfind("/")+1:]
            output_file_label2='trimmed_'+read_files[1][read_files[1].rfind("/")+1:]
            sickle_command+='pe -f '+read_files[0]+' -r '+read_files[1]+' -t sanger -o '+output_path+'/temp/'+output_file_label1\
                         +' -p '+output_path+'/temp/'+output_file_label2+' -s '+output_path+'/temp/singletons.fastq -l 100'
            check_output=output_path+'/temp/'+output_file_label1
            print(sickle_command)
            log('Read QC:\n'+sickle_command+'\n',output_path+'/log.txt')
    else:
        if len(read_files)!=1:
            print('Single-end reads requires 1 read file entered.')
            return False
        else:
            output_file_label='trimmed_'+read_files[0][read_files[0].rfind("/")+1:]
            sickle_command+='se -f '+read_files[0]+' -t sanger -o '+output_path+'/temp/'+output_file_label+' -l 100'
            check_output=output_path+'/temp/'+output_file_label
            print(sickle_command)
            log('Read QC:\n'+sickle_command+'\n',output_path+'/log.txt')
    os.system(sickle_command)
    if os.path.getsize(check_output)>0:
        return True
    else:
        return False

#RUN SPAdes
def run_spades(output_path,r_type='pe',*read_files):
    spades_command='spades.py'
    if r_type=='pe':

        if len(read_files)!=2:
            print('Paired-end reads requires 2 read files entered.')
            return False
        else:
            if check_for_file(output_path+'/temp/singletons.fastq')==True:
                spades_command+=' -k 33,55,77,99,127 -t '+str(use_x_threads)+' --only-assembler -1 '+read_files[0]+' -2 '+read_files[1]+' -s '+output_path+'/temp/singletons.fastq -o '+output_path+'/temp/assembly/'
            else:
                spades_command+=' -k 33,55,77,99,127 -t '+str(use_x_threads)+' --only-assembler -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/assembly/'
            print(spades_command)
            log('Assembly:\n'+spades_command+'\n',output_path+'/log.txt')
    else:
        if len(read_files)!=1:
            print('Single-end reads requires 1 read file entered.')
            return False
        else:
            spades_command+=' -k 33,55,77,99,127 -t '+str(use_x_threads)+' --only-assembler -s '+read_files[0]+' -o '+output_path+'/temp/assembly/'
            print(spades_command)
            log('Assembly:\n'+spades_command+'\n',output_path+'/log.txt')
    os.system(spades_command)
    if check_for_file(output_path+'/temp/assembly/scaffolds.fasta')==False:
        if check_for_file(output_path+'/temp/assembly/K127/scaffolds.fasta')==False:
            return False
        else:
            #pick K127 and copy
            shutil.copy(output_path+'/temp/assembly/K127/scaffolds.fasta',output_path+'/temp/assembly/scaffolds.fasta')
            log('Error with SPAdes scaffold selection. K127 selected.\n',output_path+'/log.txt')
            return True
    else:
        return True

#FILTER OPTIONS m AND M
# filter contigs by size
def filter_contigs_by_size(output_path,contig_file,condition,size):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    outputfile=open(contig_file,'w')
    counter=0
    for i in sequences:
        if condition=='>':
            if len(i.seq)>=int(size):
                outputfile.write('>'+i.description+'\n')
                outputfile.write(str(i.seq)+'\n')
                counter+=1
        else:
            if len(i.seq)<=int(size):
                outputfile.write('>'+i.description+'\n')
                outputfile.write(str(i.seq)+'\n')
                counter+=1
    outputfile.close()
    log('Filter by size... '+str(condition)+'='+str(size)+': '+str(counter),output_path+'/log.txt')
    del sequences[:]
    return contig_file

#CALCULATE COVERAGE USING BBMAP
def calculate_coverage(output_path,contig_file,type_reads,*read_files):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    if type_reads=='pe':       
        read1=read_files[0]
        read2=read_files[1]
        log('Read files to calculate coverage: '+read1+' & '+read2,output_path+'/log.txt')
        bbwrap_command='bbwrap.sh ref='+contig_file+' in='+read1+' in2='+read2+' out='+output_path+'/temp/aln.sam.gz kfilter=22 subfilter=15 maxindel=80'
        os.system(bbwrap_command)
        print(bbwrap_command)
        if check_for_file(output_path+'/temp/aln.sam.gz')==False:
            print('Error with BBMap. Coverage cannot be calculated.')
            log('Error with BBMap. Coverage cannot be calculated.',output_path+'/log.txt')
            return False
        pileup_command='pileup.sh in='+output_path+'/temp/aln.sam.gz '+'out='+output_path+'/temp/covstats.txt twocolumn=t'
        os.system(pileup_command)
        print(pileup_command)
        if check_for_file(output_path+'/temp/covstats.txt')==False:
            print('Error with pileup. Coverage cannot be calculated.')
            log('Error with pileup. Coverage cannot be calculated.',output_path+'/log.txt')
            return False            
    else:
        log('Read file to calculate coverage: '+read_files[0],output_path+'/log.txt')
        bbwrap_command='bbwrap.sh ref='+contig_file+' in='+read_files[0]+' out='+output_path+'/temp/aln.sam.gz kfilter=22 subfilter=15 maxindel=80'
        os.system(bbwrap_command)
        print(bbwrap_command)
        if check_for_file(output_path+'/temp/aln.sam.gz')==False:
            print('Error with BBMap. Coverage cannot be calculated.')
            log('Error with BBMap. Coverage cannot be calculated.',output_path+'/log.txt')
            return False
        pileup_command='pileup.sh in='+output_path+'/temp/aln.sam.gz '+'out='+output_path+'/temp/covstats.txt twocolumn=t'
        os.system(pileup_command)
        print(pileup_command)
        if check_for_file(output_path+'/temp/covstats.txt')==False:
            print('Error with pileup. Coverage cannot be calculated.')
            log('Error with pileup. Coverage cannot be calculated.',output_path+'/log.txt')            
            return False
    return True

#CALCULATE COVERAGE STATS TO REPORT TO USER
def calc_coverage_stats(output_path,contig_file,type_reads,reads):
    if type_reads=='pe':
        if calculate_coverage(output_path,contig_file,type_reads,reads[0],reads[1])==True:
            shutil.copyfile(output_path+'/temp/covstats.txt',output_path+'/covstats.txt')
            log('Coverage statistics written to: '+output_path+'/covstats.txt',output_path+'/log.txt')
    else:
        if calculate_coverage(output_path,contig_file,type_reads,reads[0])==True:
            shutil.copyfile(output_path+'/temp/covstats.txt',output_path+'/covstats.txt')
            log('Coverage statistics written to : '+output_path+'/covstats.txt',output_path+'/log.txt')
    return True


#FILTER OPTION C
# filter contigs by coverage
def filter_by_coverage(output_path,contig_file,val,type_reads,*read_files):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    log('Filter by coverage... >='+str(val)+': '+str(len(sequences)),output_path+'/log.txt')
    del sequences[:]
    if type_reads=='pe':
        returnFlag=calculate_coverage(output_path,contig_file,type_reads,read_files[0],read_files[1])
    else:
        returnFlag=calculate_coverage(output_path,contig_file,type_reads,read_files[0])
    if returnFlag==False:
        print('Error with pileup. Filter by coverage not performed.')
        log('Error with pileup. Filter by coverage not performed.',output_path+'/log.txt')
        return contig_file 

    #PARSE and write to file only those contigs with a coverage value >= val
    val=float(val)
    meet_threshold=[]
    cov_file=output_path+'/temp/covstats.txt'
    with open(cov_file) as csvfile:
	reader=csv.DictReader(csvfile,delimiter='\t')
	for row in reader:
            if float(row['Avg_fold'])>=val:
                meet_threshold.append(row)
    outputfile=open(contig_file,'w')
    for i in meet_threshold:
        outputfile.write('>'+i['#ID']+'\n')
        for j in sequences:
            if i['#ID']==j.id:
                outputfile.write(str(j.seq)+'\n')
                break
    outputfile.close()
    log('Filter by coverage... >='+str(val)+': '+str(len(meet_threshold)),output_path+'/log.txt')
    return contig_file

#FILTER OPTION COV
# filter contigs by spades cov value
def filter_by_spades_cov(output_path,contig_file,val):
    val=float(val)
    meet_threshold=[]
    no_cov_val=False
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    for i in sequences:
        try:
            covScore=i.id[i.id.index('_cov_'):]
            covScore=covScore[5:]
            covScore=float(covScore)
            if covScore>=val:
                meet_threshold.append(i)
        except:
            no_cov_val=True
    if no_cov_val==True:
        print('cov value not found in contigs. Contigs were not filtered.')
    else:
        outputfile=open(contig_file,'w')
        for i in meet_threshold:
            outputfile.write('>'+i.description+'\n')
            outputfile.write(str(i.seq)+'\n')
        outputfile.close()
        log('Filter by SPAdes cov value... >='+str(val)+': '+str(len(meet_threshold)),output_path+'/log.txt')
    del sequences[:]
    del meet_threshold[:]
    return contig_file

#FILTERS THE CONTIGS ACCORDING TO THE FILTER OPTIONS
# current options: min_size, max_size, min_coverage, min_cover_spades
# filter parameters are listed in variable filter_params separated by pipe
# parameter codes m=min_size; M=max_size; c=coverage; s=min_cover_spades
def filter_contigs(contig_file,filter_params,output_path,type_reads,*reads):
    filtered_contig_file=contig_file[:contig_file.rfind('/')+1]+'filtered_'+contig_file[contig_file.rfind('/')+1:]
    shutil.copy(contig_file,filtered_contig_file) #makes a copy of the file renamed with "filtered_"
    print('\nFiltering reads...')
    params=filter_params.split('|')
    for i in params:
        code,val=i.split('=')
        if code=='m':
            print('... Filter out contigs < '+val+' nucleotides in length.')
            filtered_contig_file=filter_contigs_by_size(output_path,filtered_contig_file,'>',val)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        elif code=='M':
            print('... Filter out contigs > '+val+' nucleotides in length.')
            filtered_contig_file=filter_contigs_by_size(output_path,filtered_contig_file,'<',val)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        elif code=='c':
            print('... Filter out contigs with coverage < '+val+'x.')
            filtered_contig_file=filter_by_coverage(output_path,filtered_contig_file,val,type_reads,reads[0][0],reads[0][1])
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        elif code=='cov':
            print('... Filter out contigs with a SPAdes cov value < '+val+'x.')
            filtered_contig_file=filter_by_spades_cov(output_path,filtered_contig_file,val)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
    return filtered_contig_file

#CHECKS THAT THE OUTPUT PATH EXISTS. IF IT ALREADY EXISTS, IT MAKES ANOTHER ONE.
def check_output_path(output_path):
    #check that output_path exists; if not, make it
    if not os.path.exists(output_path):         ## does not exist
        os.makedirs(output_path+'/temp')
        if not os.path.exists(output_path+'/temp'):
            print('Cannot write to output folder. Please check the output path and permissions.')
            return False,''
    else:   #ouput_path exists
        counter=1
        while os.path.exists(output_path+str(counter))==True:
            counter+=1
        print("The output path "+output_path+" already existed.\nYour output will be written to: "+output_path+str(counter)+'\n')
        output_path=output_path+str(counter)
        os.makedirs(output_path+'/temp')
        if not os.path.exists(output_path+'/temp'):
            print('Cannot write to output folder. Please check the output path and permissions.')
            return False,''
    return True,output_path


#PARSE THE BLAST AND RETURN A LIST OF DICTIONARIES
def parse_blast(filename,headers):
    x=[]
    blast_results=open(filename,'r')
    rows=csv.DictReader(blast_results,headers,delimiter=',')
    for row in rows:
        x.append(row)
    blast_results.close()
    return x

#PARSE AND PRUNE THE BLAST AND RETURN A LIST OF DICTIONARIES
def parse_N_prune_blast(file_name,headers):
    blast_results=open(file_name,'r')
    rows=csv.DictReader(blast_results,headers,delimiter=',')
    hits=[]
    for row in rows:
        hits.append(row)
    blast_results.close()

    # prune results - highest bit score wins
    best_hits=[]
    for i in hits:
        unique=True
        for j in best_hits:
            if i['qseqid']==j['qseqid']:
                if int(i['qstart'])<=int(j['qstart'])<=int(i['qend']) or int(j['qstart'])<=int(i['qstart'])<=int(j['qend']):
                    if float(i['bitscore'])>=float(j['bitscore']):
                        j['bitscore']=0
                    else:
                        unique=False
        if unique==True:
            best_hits.append(i)
        best_hits=[j for j in best_hits if not j['bitscore']==0]
    del hits[:]
    return best_hits

def filter_blast_by_min_length(blast_results,min_length):
    best_hits=[]
    for i in blast_results:
        length=abs(int(i['qstart'])-int(i['qend']))
        if length>=min_length:
            best_hits.append(i)
    return best_hits


## Run rRNA prediction
## Writes out the results to a csv file in output_path/sample_name/temp/sample_name_rRNA_predictions.csv
## output file needs to be written with sequence in it
def run_rRNA(contig_file,output_path):
    if check_for_file(database_path+'/5rnadatabase.nin')==False or check_for_file(database_path+'/16rnadatabase.nin')==False or check_for_file(database_path+'/23rnadatabase.nin')==False:
        print('DATABASE ERROR: rRNA databases could not be located. rRNA prediction not performed.\n')
        log('DATABASE ERROR: rRNA databases could not be located. rRNA prediction not performed.\n',output_path+'/log.txt')
        return False
    
    headers=['qseqid','sseqid','qstart','qend','bitscore'] #ouput format "10 qseqid sseqid qstart qend bitscore"
    #print('5S')
    output_file=output_path+'/temp/5rnadatabase.blastn'
    blastn_cline=NcbiblastnCommandline(query=contig_file, db=database_path+'/5rnadatabase', num_threads=use_x_threads, outfmt='"10 qseqid sseqid qstart qend bitscore"', out=output_file)
    stdout, stderr = blastn_cline()
    results5S=parse_N_prune_blast(output_file,headers)
    results5S=filter_blast_by_min_length(results5S,65)

    #print('16S')
    output_file=output_path+'/temp/16rnadatabase.blastn'
    blastn_cline=NcbiblastnCommandline(query=contig_file, db=database_path+'/16rnadatabase', num_threads=use_x_threads, outfmt='"10 qseqid sseqid qstart qend bitscore"', out=output_file)
    stdout, stderr = blastn_cline()
    results16S=parse_N_prune_blast(output_file,headers)
    results16S=filter_blast_by_min_length(results16S,900)

    #print('23S')
    output_file=output_path+'/temp/23rnadatabase.blastn'
    blastn_cline=NcbiblastnCommandline(query=contig_file, db=database_path+'/23rnadatabase', num_threads=use_x_threads, outfmt='"10 qseqid sseqid qstart qend bitscore"', out=output_file)
    stdout, stderr = blastn_cline()
    results23S=parse_N_prune_blast(output_file,headers)
    results23S=filter_blast_by_min_length(results23S,1800)

    # write out results to csv file
    output_file=output_path+'/temp/rRNA_predictions.csv'
    if len(results5S)>0:
        keys=results5S[0].keys()
    elif len(results16S)>0:
        keys=results16S[0].keys()
    elif len(results23S)>0:
        keys=results23S[0].keys()
    else:
        log('No rRNAs were found in the assembly\n',output_path+'/log.txt')
        print('No rRNAs were found in the assembly\n')
        return False
    with open(output_file,'w') as output_file:
        dict_writer=csv.DictWriter(output_file,keys)
        dict_writer.writeheader()
        dict_writer.writerows(results5S)
        dict_writer.writerows(results16S)
        dict_writer.writerows(results23S)
    del results5S[:]
    del results16S[:]
    del results23S[:]

    #write out rRNA sequences to file
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    rrna_sequences=[]

    list_rrnas=[]
    with open(output_path+'/temp/rRNA_predictions.csv','rb') as infile:
	reader=csv.DictReader(infile)
	for row in reader:
		list_rrnas.append(row)
    
    for i in list_rrnas:
        for k in sequences:
            if i['qseqid']==k.description:
                start=int(i['qstart'])
                stop=int(i['qend'])

                if start>stop:
                    #it's reverse complement
                    subseq=k.seq[stop-1:start].reverse_complement()
                else:
                    subseq=k.seq[start-1:stop]
                molecule=i['sseqid'][i['sseqid'].find('|')+1:i['sseqid'].find('|',i['sseqid'].find('|')+1)]
                rrna_label='>'+i['qseqid']+'_'+molecule+'_'+str(start)+'_'+str(stop)
                rrna_sequences.append([rrna_label,subseq])
                break
    rrna_output=open(output_path+'/predicted_rRNAs.fasta','w')
    for i in rrna_sequences:
        rrna_output.write(str(i[0])+'\n')
        rrna_output.write(str(i[1])+'\n')
    rrna_output.close()
    del rrna_sequences[:]
    return True

#check orfs to return only those that do not overlap with RNA predictions
def check_orfs(output_path):
    orf_seq_file=output_path+'/temp/predicted_orfs.fasta'
    rRNA_file=output_path+'/temp/rRNA_predictions.csv'
    tRNA_file=output_path+'/temp/tRNAscan_locations.csv'

    #make list of RNA locations
    list_rrnas=[]
    if check_for_file(rRNA_file)==True:
        with open(rRNA_file,'rb') as infile:
            reader=csv.DictReader(infile)
            for row in reader:
                    list_rrnas.append(row)
    #need to add tRNAs
    list_trnas=[]
    if check_for_file(tRNA_file)==True:
        with open(tRNA_file,'rb') as infile:
            reader=csv.DictReader(infile)
            for row in reader:
                    list_trnas.append(row)		

    orf_sequences=list(SeqIO.parse(orf_seq_file,'fasta'))
    final_orfs=[]
    for i in orf_sequences:
        orf_header=i.id
        start=int(orf_header[orf_header.rfind('_',1,orf_header.rfind('_')-1)+1:orf_header.rfind('_')])
        end=int(orf_header[orf_header.rfind('_')+1:])
        contig_name=orf_header[1:orf_header.find('_orf')]
        orf_start=min(start,end)
        orf_end=max(start,end)

        conflict=False
        for j in list_rrnas:
            if contig_name==j['qseqid']:
                rna_start=min(int(j['qstart']),int(j['qend']))
                rna_end=max(int(j['qstart']),int(j['qend']))

                if rna_start <= orf_start <=rna_end or rna_start <= orf_end <=rna_end:
                    conflict=True
        if conflict==False:
            for k in list_trnas: #'tRNA,contig,start,end,tRNAscan-SE_score')
                if contig_name==k['contig']:
                    trna_start=min(int(k['start']),int(k['end']))
                    trna_end=max(int(k['start']),int(k['end']))
                    if trna_start <= orf_start <=trna_end or trna_start <= orf_end <=trna_end:
                        conflict=True                    
        if conflict==False:
            final_orfs.append(i)
    SeqIO.write(final_orfs,output_path+'/predicted_orfs_nt.fasta','fasta')

    orf_output=open(output_path+'/predicted_orfs_aa.fasta','w')
    for i in final_orfs:
        orf_output.write('>'+str(i.id)+'\n')
        orf_output.write(str(i.seq.translate())+'\n')
    orf_output.close()

    #log results
    log('# final orfs: '+str(len(final_orfs))+'\n',output_path+'/log.txt')
    print('Final # of ORFs: '+str(len(final_orfs))+'\n')

    del final_orfs[:]
    del orf_sequences[:]
    del list_rrnas[:]
    return True

#predicts the functions of ORFs based upon user supplied databases.
def predict_function(output_path,database_path,database_name):
    output_file=output_path+'/temp/func_orfs.tblastn'
    query_file=output_path+'/predicted_orfs_aa.fasta'

    #get names of all orfs
    orfs=[]
    seqs=list(SeqIO.parse(query_file,'fasta'))
    for i in seqs:
        orfs.append(i.id)
    del seqs[:]

    print('Predicting ORF functions based on homology.\nComparing to database: '+str(database_name[0])+'.\n')
    tblastn_cline=NcbitblastnCommandline(query=query_file, db=database_path+'/'+str(database_name[0]), num_threads=use_x_threads, max_target_seqs=1, outfmt='"10 qseqid sseqid qcovs pident length evalue bitscore"', out=output_file)
    stdout, stderr = tblastn_cline()

    headers=['qseqid','sseqid','qcovs','pident','length','evalue','bitscore']
    b=parse_blast(output_file,headers)

    annotation=[]
    #if meets threshold, add to the annotation list
    for i in b:
        if float(i['qcovs'])>homol_qcov and float(i['pident'])>homol_pident and float(i['bitscore'])>homol_bitscore:
            i['db']=str(database_name[0])
            annotation.append(i)
    annotation_includes=[i['qseqid'] for i in b]

    #search through other databases
    otherdb=database_name[1:]
    for i in otherdb:
        print('Predicting ORF functions based on homology.\nComparing to database: '+i+'.\n')
        tblastn_cline=NcbitblastnCommandline(query=query_file, db=database_path+'/'+i, num_threads=use_x_threads, max_target_seqs=1, outfmt='"10 qseqid sseqid qcovs pident length evalue bitscore"', out=output_file)
        stdout, stderr = tblastn_cline()
        b=parse_blast(output_file,headers)

        #determine if better or not in annotation -- no -- i don't think it'll be in the same order? is this doing what i think it is doing?
        for j in b:
            try:
                index=annotation_includes.index(j['qseqid'])
            except:
                index=-1
            
            if index==-1:
                j['db']=i
                annotation.append(j)
            else:
                #is it better than the annotation selected now?
                if float(j['qcovs'])>annotation[index]['qcovs'] and float(j['pident'])>annotation[index]['pident'] and float(j['bitscore'])>annotation[index]['bitscore']:
                    j['db']=1
                    del annotation[index]
                    annotation.append(j)
        annotation_includes=[j['qseqid'] for j in b]
    del b[:]

    #write out the annotations
    outputf=open(output_path+'/annotations.csv','w')
    #output format contig#,orfID,start,end,gene_name,synonym,COG,product,homolog_info
    outputf.write('contig,orfID,start,end,gene_name,synonym,COG,product,homolog_info\n')
    for i in orfs:
        orf_header=i
        start=int(orf_header[orf_header.rfind('_',1,orf_header.rfind('_')-1)+1:orf_header.rfind('_')])
        end=int(orf_header[orf_header.rfind('_')+1:])
        contig_name=orf_header[:orf_header.find('_orf')]
        orf_name=orf_header[len(contig_name)+1:orf_header.rfind('_',1,orf_header.rfind('_')-1)]

        outputf.write(contig_name+','+orf_name+','+str(start)+','+str(end)+',')
        found=False
        for j in annotation:
            if j['qseqid']==orf_header and found==False:
                #output the info
                homolog=j['sseqid']
                gene_name=homolog[:homolog.find('|')]
                synonym=homolog[len(gene_name)+1:homolog.find('|',len(gene_name)+1)]
                cog=homolog[len(gene_name)+1+len(synonym)+1:homolog.find('|',len(gene_name)+1+len(synonym)+1)]
                product=homolog[len(gene_name)+1+len(synonym)+1+len(cog)+1:homolog.find('|',len(gene_name)+1+len(synonym)+1+len(cog)+1)]
                homolog_info=homolog[len(gene_name)+len(synonym)+len(cog)+len(product)+4:]
                outputf.write(gene_name+','+synonym+','+cog+','+product+','+homolog_info+'\n')
                found=True
        if found==False:  #no homology detected, output as hypothetical protein
            outputf.write('-,-,-,PRED:HYPOTHETICAL PROTEIN,-\n')

    outputf.close()
    del annotation[:]
    del annotation_includes[:]
    return True

#run tRNA
def run_tRNA(contig_file,output_path):
    #trnascan_command='tRNAscan-SE -B '+contig_file+' -o '+output_path+'/tRNAscan-SE_output.txt'
    trnascan_command='tRNAscan-SE -o '+output_path+'/tRNAscan-SE_output.txt -B '+contig_file
    print(trnascan_command)
    os.system(trnascan_command)

    #write out temp file of locations
    with open(output_path+'/tRNAscan-SE_output.txt') as f:
        lines=f.read().splitlines()
    if len(lines)<4:
        print('No tRNAs found.')
        return 0
    lines=lines[3:] #strip header
    fOut=[] #new list to hold formatted output
    fOut.append('tRNA,contig,start,end,tRNAscan-SE_score')
    for l in lines:
        if l!='':
            l=l.replace(' ','') #kill whitespaces
            l=l.split('\t')
            l=','.join((l[4],l[0],l[2],l[3],l[8]))
            fOut.append(l)
    out_file=open(output_path+'/temp/tRNAscan_locations.csv','w')
    out_file.write('\n'.join(fOut))
    out_file.close()
    return len(lines)

#rename contigs to consistent numbering after processing
def rename_contigs(contig_file,output_path):
    contigs=[]
    contigs=list(SeqIO.parse(contig_file,'fasta'))
    outputf=open(output_path+'/final_contigs.fasta','w')
    counter=1
    for i in contigs:
        i.name=''
        i.description=''
        i.id='contig_'+str(counter)
        counter=counter+1
    SeqIO.write(contigs,output_path+'/final_contigs.fasta','fasta')
    del contigs[:]
    return output_path+'/final_contigs.fasta'

#get number of sequences in a FASTA file
def get_num_seqs(seq_file):
    contigs=list(SeqIO.parse(seq_file,'fasta'))
    num_contigs=len(contigs)
    return num_contigs


#RUN PIPELINE
def run_peasant(assembly,contig_file,type_reads,reads,filter_params,output_path,database_name):
    if assembly!='provided':
        if qc_and_assembly(output_path,assembly,type_reads,reads)==False:
            print('\nError in QC and assembly. peasant did not run.')
            log('Error in QC and assembly. peasant did not run.\n',output_path+'/log.txt')
            return False
        else:
            contig_file=output_path+'/temp/assembly/scaffolds.fasta'
            if type_reads=='pe':
                log('Assembly produced by '+assembly+'.\nRead file(s):\n\t'+reads[0]+'\n\t'+reads[1]+'\n',output_path+'/log.txt')
            else:
                log('Assembly produced by '+assembly+'.\nRead file(s):\n\t'+reads+'\n',output_path+'/log.txt')
    else: #copy assembly file into output_path/temp/
        log('Assembly file provided: '+contig_file+'\n',output_path+'/log.txt')
        temp_output_file=output_path+'/temp/'+contig_file[contig_file.rfind("/")+1:]
        shutil.copyfile(contig_file,temp_output_file)
        contig_file=temp_output_file
    num_contigs_from_assembly=get_num_seqs(contig_file)
    log('# contigs in provided assembly: '+str(num_contigs_from_assembly)+'\n',output_path+'/log.txt')

    #RUN FILTERS
    if num_contigs_from_assembly==0:
        print('Error with assembly. No contigs were created.')
        log('Error with assembly. No contigs were created.\n',output_path+'/log.txt')
        return False
    contig_file=qc_contig_file(contig_file) #reformat contigs file
    if contig_file==False:
        print('\nError with QC of contig file. Analysis cannot continue.')
        log('Error with QC of contig file. Analysis cannot continue.\n',output_path+'/log.txt')
        return False
    print('Assembly produced '+str(num_contigs_from_assembly)+' contigs')
    
    if filter_params!='':
        # run filters
        contig_file=filter_contigs(contig_file,filter_params,output_path,type_reads,reads)
        num_contigs_after_filtering=get_num_seqs(contig_file)
        print('# contigs after filtration: '+str(num_contigs_after_filtering))
        log('# contigs after filtration: '+str(num_contigs_after_filtering)+' contigs',output_path+'/log.txt')
        
        if num_contigs_from_assembly==num_contigs_after_filtering:
            print('Filters did not remove any contigs.')
            log('Filters did not remove any contigs.\n',output_path+'/log.txt')           
        else:
            print('Filters removed '+str(num_contigs_from_assembly-num_contigs_after_filtering)+' contigs.')
            log('Filters removed '+str(num_contigs_from_assembly-num_contigs_after_filtering)+' contigs.\n',output_path+'/log.txt')
    else:
        print('No filters applied.')
        log('No filters applied.\n',output_path+'/log.txt')

    log('N50 of filtered assembly: '+str(N50(contig_file))+'\n',output_path+'/log.txt')
    log('Number of nucleotides in filtered assembly: '+str(assembly_size(contig_file))+'\n',output_path+'/log.txt')
    log('GC content of filtered assembly: '+str(calc_GC_content(contig_file))+'\n',output_path+'/log.txt')
    calc_coverage_stats(output_path,contig_file,type_reads,reads)
    

    ##RUN GLIMMER AND PREDICTIONS
    contig_file=rename_contigs(contig_file,output_path)
    orf_file=run_glimmer(contig_file,output_path)
    num_predicted_orfs=get_num_seqs(orf_file)
    print('# ORFs predicted by GLIMMER: '+str(num_predicted_orfs))
    log('# ORFs predicted by GLIMMER: '+str(num_predicted_orfs),output_path+'/log.txt')
    if num_predicted_orfs==0:
        print('No ORFs were predicted')
        log('No ORFs were predicted.\n',output_path+'/log.txt')
        return False
    if run_rRNA(contig_file,output_path)==True:
        num_predicted_rRNAs=get_num_seqs(output_path+'/predicted_rRNAs.fasta')
        print('# rRNAs predicted: '+str(num_predicted_rRNAs))
        log('# rRNAs predicted: '+str(num_predicted_rRNAs),output_path+'/log.txt')        
        
    num_tRNAs=run_tRNA(contig_file,output_path)
    log('# tRNAs predicted: '+str(num_tRNAs)+'\n',output_path+'/log.txt')
    check_orfs(output_path)
    log('Final ORF sequence files created.\n',output_path+'/log.txt')
    predict_function(output_path,database_path,database_name)
    return True




############# THIS IS ACTUALLY DOING ALL OF THE CHECKS AND RUNNING IT #######################

ts=time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
message2log='Running peasant.\n\n\n'+st+'\n'
command='python'
for i in sys.argv:
    command+=' '+i
message2log=message2log+command


#SET VARIABLES
if args.assembler is None:
    assembly='provided'
    contig_file=args.assembled_contigs
    type_reads=None
    reads=None
else:
    contig_file=''
    assembly=args.assembler
    if args.single_reads is not None:
        type_reads='se'
        reads=args.single_reads
    else:
        type_reads='pe'
        reads=args.paired_end_reads
output_path=args.output_path

if args.num_threads is None:
    use_x_threads=1
else:
    use_x_threads=args.num_threads
if args.num_threads>multiprocessing.cpu_count():
    use_x_threads=1

filter_params=''
if args.min_contig_size is not None:
    filter_params+='m='+str(args.min_contig_size)
if args.max_contig_size is not None:
    if len(filter_params)>0:
        filter_params+='|'
    filter_params+='M='+str(args.max_contig_size)
if args.min_coverage is not None:
    if len(filter_params)>0:
        filter_params+='|'
    filter_params+='c='+str(args.min_coverage)
if args.min_SPAdes_cov is not None:
    if len(filter_params)>0:
        filter_params+='|'
    filter_params+='cov='+str(args.min_SPAdes_cov)
database_name=args.genus
if args.qcov is None:
    homol_qcov=70.0
if args.pident is None:
    homol_pident=70.0
if args.bitscore is None:
    homol_bitscore=50.0
message2log=message2log+'\nHomology thresholds:\n\tquery coverage='+str(homol_qcov)+'\n\t%ID='+str(homol_pident)+'\n\tbitscore='+str(homol_bitscore)+'\n'


#CHECK PATH, FILES, and PATHS PROVIDED BY THE USER
x=check_installs() #check that the PATH variables are set correctly
if x==False:
    print('Path variables must be set prior to execution.\nCheck PATH and try again.')
else:
    #check reads file (if starting from raw data); check assembly if starting from there
    if assembly=='provided':
        fileQC,errmsg=check_files(contig_file)
    else:
        if type_reads=='se':
            fileQC,errmsg=check_files(reads)
        else:
            fileQC,errmsg=check_files(reads[0],reads[1])

    if fileQC==False:   #handles errors due to reads or assembly
        print(errmsg)
    else:               #reads or assembly a-okay
        #check DB locations
        dbfilesokay=True
        for i in database_name:
            fileQC,errmsg=check_files(database_path+'/'+i+'.nhr')
            if fileQC==False:
                #check to see if it is a db alias:
                fileQC,errmsg_alias=check_files(database_path+'/'+i+'.nal')
                if fileQC==False:
                    print('DATABASE ERROR: '+errmsg)
                    dbfilesokay=False
        if dbfilesokay:
            x,output_path=check_output_path(output_path)    #makes output path
            log(message2log,output_path+'/log.txt')
            if x==True:
                runOK=run_peasant(assembly,contig_file,type_reads,reads,filter_params,output_path,database_name)
                if runOK==True:
                    print('Your job has finished running.\n')
                    log('Your job has finished running.\n',output_path+'/log.txt')
                    ts=time.time()
                    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
                    log('Job completed at: '+st+'\n',output_path+'/log.txt')
                    print('Log written: '+output_path+'/log.txt\n')
                else:
                    print('Error occurred. Check parameters and files and try again.')
                    log('Error occurred. Check parameters and files and try again.',output_path+'/log.txt')
            else:
                print('Could not locate/create the output directory.')
