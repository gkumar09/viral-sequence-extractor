#!/usr/bin/env python

####################
# running command: ./pipeline.py seq_data1_ncbi/genomes_PP seqData
#
#####################

import os
import sys
import subprocess
import glob
import shutil

subprocess.call('pip install pysam==0.8.4', shell=True)
subprocess.call('pip install pysamstats', shell=True) 

import pysam
import pysamstats

print ('The current working directory is: ' + os.getcwd())

# bowtie index building
fullgenome= subprocess.call ("cat " + sys.argv[1] +  "/* >> fullgenome.fa", shell= True)
subprocess.call("bowtie2-2.2.6/bowtie2-build fullgenome.fa fullgenome", shell=True)

os.mkdir('fullgenome')
for files in glob.glob("*.bt2"):
    shutil.move(files, "fullgenome")

# reading manifest file and running bowtie2
fil= open("manifestfile.csv")
files=fil.readlines()[1:]
for reads in files: 
    pereads= reads.split(',')
  
    # run bowtie
    print ('The command used: ' + 'bowtie2-2.2.6/bowtie2 fullgenome.fa -1' + ' ' + sys.argv[2]+ '/' + pereads[2] + " " + "-2 " + sys.argv[2] + '/' + pereads[3] + ' '+ '-S' + ' '+ 'bowtie_output_'+pereads[0]+'.sam')
    print (' ')
    subprocess.call('bowtie2-2.2.6/bowtie2 -x fullgenome/fullgenome -1' + ' ' + sys.argv[2]+ '/' + pereads[2] + " " + "-2 " + sys.argv[2] + '/' + pereads[3] + ' '+ '-S' + ' '+ 'bowtie_output_'+pereads[0]+'.sam', shell=True)
  
    # get consensus file
    print ('The command used: ' + 'samtools-1.3/samtools view -bS' + ' ' + 'bowtie_output_'+pereads[0]+'.sam ' + '>>' + ' bowtie_output_'+pereads[0]+'.bam') 
    subprocess.call('samtools-1.3/samtools view -bS' + ' ' + 'bowtie_output_'+pereads[0]+'.sam ' + '>>' + ' bowtie_output_'+pereads[0]+'.bam', shell=True)
 
    #print ('The command used: ' + 'samtools-1.3/samtools sort' + ' ' + 'bowtie_output_'+pereads[0]+'.bam ' + '-o' + ' bowtie_output_'+pereads[0]+'_sorted.bam')
    subprocess.call('samtools-1.3/samtools sort' + ' ' + 'bowtie_output_'+pereads[0]+'.bam ' + '-o' + ' bowtie_output_'+pereads[0]+'_sorted.bam', shell=True) 

    #print ('The command used: ' + 'samtools-1.3/samtools mpileup -uf' + ' ' + 'fullgenome.fa' + ' bowtie_output_'+pereads[0]+'_sorted.bam ' + '| ' + 'bcftools-1.3/bcftools call -mv -Oz -o ' + 'raw_calls_'+pereads[0]+'.vcf.gz')
    subprocess.call('samtools-1.3/samtools mpileup -uf' + ' ' + 'fullgenome.fa' + ' bowtie_output_'+pereads[0]+'_sorted.bam ' + '| ' + 'bcftools-1.3/bcftools call -mv -Oz -o ' + 'raw_calls_'+pereads[0]+'.vcf.gz', shell=True)

    #print ('The command used: ' + 'htslib-1.3/tabix' + ' ' + 'raw_calls_'+pereads[0]+'.vcf.gz')  
    subprocess.call('htslib-1.3/tabix' + ' ' + 'raw_calls_'+pereads[0]+'.vcf.gz', shell=True)

    #print ('The command used: ' + 'gunzip' + ' ' + 'raw_calls_'+pereads[0]+'.vcf.gz')                                                                                   
    subprocess.call('gunzip' + ' ' + 'raw_calls_'+pereads[0]+'.vcf.gz', shell=True)
    
    subprocess.call('grep' + ' ' + '-A100000' + ' ' + '"#CHROM"' + ' raw_calls_'+pereads[0]+'.vcf' + '|' + 'cut -f1,2,4,5>>' + pereads[0]+'_consensus_results.txt', shell=True)

    # index the .bam file
    print ('The command used: ' + 'samtools-1.3/samtools index' + ' ' + 'bowtie_output_'+pereads[0]+'_sorted.bam')
    subprocess.call('samtools-1.3/samtools index' + ' ' + 'bowtie_output_'+pereads[0]+'_sorted.bam', shell=True)

    print ('The command used: ' + 'pysamstats -f ' + ' fullgenome.fa --type variation_strand' + ' bowtie_output_'+pereads[0]+'_sorted.bam ' + '--format=csv >> ' + pereads[0]+'_nucleotide')
    subprocess.call('pysamstats -f ' + ' fullgenome.fa --type variation_strand' + ' bowtie_output_'+pereads[0]+'_sorted.bam ' + '--format=csv >> ' + pereads[0]+'_nucleotide', shell=True)    

    print ('The command used: ' + 'cat ' + pereads[0]+'_nucleotide' + '| sort | uniq >> ' + pereads[0]+'_nucleotidepercentage')
    subprocess.call('cat ' + pereads[0]+'_nucleotide' + '| sort | uniq >> ' + pereads[0]+'_nucleotidepercentage', shell=True)

    consen_file= open(pereads[0]+'_consensus_results.txt')
    consensus_file= consen_file.readlines()[1:]
   
    nucleo_file= open(pereads[0]+'_nucleotidepercentage')
    nucleotide_file= nucleo_file.readlines()[:-1]
   
    out_file= open(pereads[0]+'_final_percentage_nucleotides', "w")
    
    def consensus_nucleotide():
        # processing first file
        filt_ff=[]
        for filt1 in consensus_file:
            filt1_tab= filt1.replace('\t', ',')
            filt1_nline= filt1_tab.split('\n')
            filt_ff.append(filt1_nline)
    
        filt2_ff=[]
        for filt2 in range(0, len(filt_ff)):
            filt2_sep= filt_ff[filt2][0]
            filt2_com= filt2_sep.split(',')
            filt2_join= ','.join(filt2_com[:2])
            filt2_ff.append(filt2_join) 
    
        # processing second file
        filt_sf=[]
        for filt1_2 in nucleotide_file:
            filt1_2_nline= filt1_2.split('\n')
            filt_sf.append(filt1_2_nline)
    
        for filt2_2 in range(0, len(filt_sf)):
            filt2_2_out= filt_sf[filt2_2][0]
            filt2_2_com= filt2_2_out.split(',')
            filt2_2_join= ','.join(filt2_2_com[0:2])
            for nucl in range(len(filt2_ff)):
                if filt2_2_join == filt2_ff[nucl]:
                    nucl_per= filt2_2_out
                    nucl_per_com= nucl_per.split(',')
                    out_file.write("\n")
                    out_file.write(str('At the genome and position:') + str(filt2_ff[nucl]) + '\n')
                    out_file.write(str('The percentage of A is: ') + str(round(float(float(nucl_per_com[33])/float(nucl_per_com[3]))*100, 2))+ '\n')
                    out_file.write(str('The percentage of C is: ') + str(round(float(float(nucl_per_com[39])/float(nucl_per_com[3]))*100, 2))+ '\n')
                    out_file.write(str('The percentage of T is: ') + str(round(float(float(nucl_per_com[45])/float(nucl_per_com[3]))*100, 2))+ '\n')
                    out_file.write(str('The percentage of G is: ') + str(round(float(float(nucl_per_com[51])/float(nucl_per_com[3]))*100, 2))+ '\n')
                    out_file.write("\n")
       
        out_file.close()
        
  
    consensus_nucleotide()


    subprocess.call('cat ' + pereads[0]+'_consensus_results.txt' + ' ' + pereads[0]+'_final_percentage_nucleotides' + ' >> ' + pereads[0]+'_consensus_nucleotide_results.txt', shell=True)


subprocess.call('rm *.bam', shell=True)
subprocess.call('rm *_sorted.bam.bai', shell=True) 
subprocess.call('rm *_consensus_results.txt', shell=True) 
subprocess.call('rm *_final_percentage_nucleotides', shell=True)
subprocess.call('rm *_nucleotidepercentage', shell=True) 
subprocess.call('rm *_nucleotide', shell=True) 
subprocess.call('rm raw_calls*', shell=True) 
subprocess.call('rm *.fa.fai', shell=True)
