#!/usr/bin/python
from oss import oss_api
from oss import oss_util
from multiprocessing import Pool
import getopt
import sys
import os
import time
import glob
try:
	Host = os.environ['ALI_DIKU_OSS_HOST']
except:
	Host = 'oss-cn-qingdao.aliyuncs.com'
	pass

Id   = 'QFmrMPB18qNx9KYc'
Key  = 'IuAdh4qL9noDf0UnMOO977HSgZSc0E'
try:
	oss = oss_api.OssAPI(Host,Id,Key)
except:
	print "oss inition error!"
	sys.exit()

bed=""
bucket  = ""
objname1 = ""
objname2 = ""
start_state = 0
local =False
download = False
program = 'varscan'
local_out_dir = '/expan/programs/vcfs'

#YUN VARIALBLES

picard_samsort = '/root/Downloads/picard-tools-1.119/SortSam.jar'
picard_deduption= '/root/Downloads/picard-tools-1.119/MarkDuplicates.jar'
picard_add_rg= '/root/Downloads/picard-tools-1.119/AddOrReplaceReadGroups.jar'
picard_create_dict= '/root/Downloads/picard-tools-1.119/CreateSequenceDictionary.jar'
picard_reordersam = '/expan/programs/fastq2vcf/picard-tools-1.119/ReorderSam.jar'
varscan = '/root/Downloads/VarScan.v2.3.9.jar'
gatk = '/root/Downloads/GenomeAnalysisTK.jar'
memory = '40g'
threads = 15
ref = '/home/hg19/hg19.fa'
tmp_dir = '/home/admin/tmp'
out_dir= 'jingyun-output'

########some thing about options##########
usage ='''
	python run_job.py [options]
		-h  --help       read help file
		-b  --bucket     bucket name better in qingdao
		-f  --file1      fastq file1
		-g  --file2      fastq file2
		-t  --threads    num of threads
		-s  --state      from where to start the program
		-l  --local      run_local version
		-d  --download   just down data specifized
		-p  --program    gatk or varscan(default)
	   	-m  --memory     defaut  40g
		-r  --bed		 target bed file
		-a  --anno	     sample_no
		-o  --out_dir    out put dir on oss
		'''
try:
	options,args = getopt.getopt(sys.argv[1:],"hb:f:g:t:s:ldp:m:r:a:o:",["help","bucket=","file1=","file2=","threads=","state=","local","download","program=","memory=","bed=","anno","out_dir="])
except:
	print usage
	sys.exit()
for name,value in options:
	if name in ('-h','--help'):
		print usage
		sys.exit()
	if name in ('-f','file1='):
		objname1 = value
	if name in ('-b','bucket='):
		bucket = value
	if name in ('-g','file2='):
		objname2 = value
	if name in ('-t','threads='):
		threads = value
	if name in ('-s','state='):
		start_state = int(value)
	if name in ('-l','--local'):
		local = True
		#LOCAL VARIABLES
		picard_samsort = '/expan/programs/fastq2vcf/picard-tools-1.119/SortSam.jar'
		picard_deduption= '/expan/programs/fastq2vcf/picard-tools-1.119/MarkDuplicates.jar'
		picard_add_rg= '/expan/programs/fastq2vcf/picard-tools-1.119/AddOrReplaceReadGroups.jar'
		picard_create_dict= '/expan/programs/fastq2vcf/picard-tools-1.119/CreateSequenceDictionary.jar'
		picard_reordersam = '/expan/programs/fastq2vcf/picard-tools-1.119/ReorderSam.jar'
		varscan = '/expan/programs/fastq2vcf/VarScan.v2.3.9.jar'
		gatk = '/expan/programs/fastq2vcf/GenomeAnalysisTK.jar'
		memory = '5g'
		threads = 5
		ref = '/expan/programs/fastq2vcf/hg19/hg19.fa'
		tmp_dir = '/expan/tmp/'
	if name in ('-d','--download'):
		download = True
	if name in ('-p','program='):
		program = value
	if name in ('-m','memory='):
		memory = value
	if name in ('-r','bed='):
		bed = value
	if name in ('-a','anno='):
		sample_no = value
	if name in ('-o','out_dir='):
		out_dir = value
###### function start ####################
def used_time(func):
	def _used_time(*args,**kargs):
		start = time.time()
		ret = func(*args,**kargs)
		end = time.time()
		print "%s used time:%d seconds" % (func.__name__,end-start)
		return ret
	return _used_time
def count_reads(bam,option,region=''):
    '''
    options is ['all','mapped','region']
    if option is region then you have to specify the region arguments like 'chr2:10,000,20,000'
    '''
    if option == 'all':
        num = os.popen('samtools view -c %s' % bam).readlines()
        num = int(num[0].strip('\n'))
        return num
    if option == 'mapped':
        num = os.popen('samtools view -c -F 4 %s ' % bam).readlines()
        num = int(num[0].strip('\n'))
        return num
    if option == 'region':
        num = os.popen('samtools view -c %s %s' % (bam,region)).readlines()
        num = int(num[0].strip('\n'))
        return num

def compute_maped_and_target_ratio_out2file(raw_bam,bed_bam):
	#compute
	all_reads_num = count_reads(raw_bam,'all')
	mapped_reads_num = count_reads(raw_bam,'mapped')
	target_reads_num = count_reads(bed_bam,'all')
	mapped_ratio = round(float(mapped_reads_num)/float(all_reads_num),2)
	target_ratio = round(float(target_reads_num)/float(all_reads_num),2)
	#store outcome to a file
	prex = raw_bam.split('.')[0]
	out = prex + '.mapped_target.ratio.txt'
	fp = open(out,'w')
	fp.write('%s\n%s\n' % (mapped_ratio,target_ratio))
	fp.close()
	return out
def split_bed(bed):
	fp = open(bed,'r')
	n=0
	j=0
	fp_out = open('0.split.bed','w')
	for line in fp.readlines():
		n= n + 1
		if n % 5000 == 0:
			j=j+1
			out = str(j) + '.split.bed'
			fp_out = open(out,'w')
			fp_out.write(line)
		else:
			fp_out.write(line)
	bed_list = glob.glob('*.split.bed')
	return bed_list
def compute_region_reads_count(bed_bam,bed):
	out = bed.split('.')[0] + '.target.txt'
	fp_out = open(out,'w')
	fp = open(bed,'r')
	lines = fp.readlines()
	lines = [x.strip('\n') for x in lines if not x == '\n']
	for line in lines:
		line = line.split('\t')
		chrome = line[0]
		start = line[1]
		end = line[2]
		gene = line[3]
		query_stry = chrome + ':' + start + '-' + end
		reads_num = count_reads(bed_bam,'region',query_stry)
		line_out = '%s\t%s\t%s\t%s\t%s\n' % (chrome,start,end,gene,reads_num)
		fp_out.write(line_out)
	fp_out.close()
@used_time
def Compute_regions_reads(bed_bam,big_bed):
	beds = split_bed(big_bed)
	pools = Pool(20)
	for bed in beds:
		pools.apply_async(compute_region_reads_count,(bed_bam,bed))
	pools.close()
	pools.join()
	os.system('rm -f *.split.bed')
	targets = glob.glob('*.target.txt')
	out = bed_bam.split('.')[0] + '.Target.txt'
	for i in range(len(targets)):
		os.system('cat %s.target.txt >> %s' % (i,out))
	os.system('rm -f *.target.txt')
	return out

def get_file_from_oss(bucket,objname):
	#Get fastq1 and fastq2 from oss
	##prepare files
	filename = objname.split('/')[-1]
	try:
		res = oss_util.multi_get(oss,bucket,objname,filename,10,20)
		if res.real == 1:
			print "Load file %s successfully!" % objname
			return filename
		else:
			print "Load file %s error!" % objname
			sys.exit()
	except:
		print "Get file from oss error!\nyou can retry or contact kongdeju for detail\n"
		sys.exit()
@used_time
def get_two_file(bucket,objname1,objname2):
	gfile1 = get_file_from_oss(bucket,objname1)
	gfile2 = get_file_from_oss(bucket,objname2)	
	return gfile1,gfile2

@used_time
def get_one_file(bucket,objname1):
	gfile1 = get_file_from_oss(bucket,objname1)
	return gfile1

def uncompress(filename):
	try:
		os.system('gzip -df %s' % filename)
		filename = filename.strip('.gz')
		return filename
	except:
		print "Uncompressing error!"
		sys.exit()

def compress(filename):
	try:
		os.system('gzip -f %s' % filename)
	except:
		print "compressing error!"
		sys.exit()

@used_time
def uncompress_two_file(gfile1,gfile2):
	qfile1 = uncompress(gfile1)
	qfile2 = uncompress(gfile2)
	return qfile1,qfile2

@used_time
def uncompress_one_file(gfile1):
	qfile1 = uncompress(gfile1)
	return qfile1

@used_time
def bwa_map(file1,file2):
	tmp_name = file1.split('.')[0]
	tmp_name = '%s.sam' % tmp_name
	os.system('bwa mem -t %s -M %s %s %s > %s' % (threads,ref,file1,file2,tmp_name))
	return tmp_name
@used_time
def picard_sort(sam):
	prex = sam.rstrip('.sam')
	tmp_name1 = prex + '.bam'
	tmp_name2 = prex + '.sort.bam'
	os.system('samtools view -bS %s > %s' % (sam,tmp_name1))
	os.system('java -Xmx%s  -jar %s INPUT=%s OUTPUT=%s SORT_ORDER=coordinate TMP_DIR=%s' % (memory,picard_samsort,tmp_name1,tmp_name2,tmp_dir))
	return tmp_name2
@used_time
def picard_dedup(bam):
	prex = bam.split('.')[0]
	tmp_name = prex + '.dedup.bam'
	os.system('java -Xmx%s -jar %s INPUT=%s OUTPUT=%s REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT M=tmp.stat TMP_DIR=%s' % (memory,picard_deduption,bam,tmp_name,tmp_dir))
	return tmp_name

def format_bed(bed):
    fp = open(bed,'r')
    bed_new = bed.split('.')[0]+'.new.bed'
    fp_new = open(bed_new,'w')
    lines = fp.readlines()
    fp.close
    for line in lines :
        line = line.strip('\n\r')
        if '\t' in line:
            line = line.replace(' ','')
        else:
            line = '\t'.join(line.split())
        line = line + '\n'
        fp_new.write(line)
    fp_new.close
    print 'new bed is %s' % bed_new
    return bed_new

@used_time
def intersect(bam,bed2):
	global bed 
	bed = format_bed(bed2)
	prex = bam.split('.')[0]
	target_bam = prex + '.target.bam'
	os.system("sed -i 's/^chr//' %s" % bed)
	os.system("sed -i 's/^M/MT/' %s" % bed)
	os.system("intersectBed -abam %s -b %s -wa -u >%s" % ( bam,bed,target_bam))
	os.system("samtools index %s" % target_bam)
	return target_bam

@used_time
def mpileup(dedup_bam):
	pre = dedup_bam.split('.')[0]
	tmp_name = pre + '.pipleup'
	os.system('samtools mpileup -f %s %s > %s' %(ref,dedup_bam,tmp_name))
	return tmp_name
@used_time
def varscan_call(pileup):
	tmp_name = pileup.split('.')[0]
	tmp_name1 = tmp_name + '.snp'
	tmp_name2 = tmp_name + '.indel'
	os.system('java -Xmx%s -jar %s mpileup2snp %s  --output-vcf 1>%s' % (memory,varscan,pileup,tmp_name1))
	os.system('java -Xmx%s -jar  %s mpileup2indel %s  --output-vcf 1 > %s' % (memory,varscan,pileup,tmp_name2))
	return tmp_name1,tmp_name2

@used_time
def snp_indel_merge(snp,indel):
	prex = snp.split('.')[0]
	snp_gz = snp + '.gz'
	indel_gz = indel + '.gz'
	con_vcf = prex + '.all.vcf'
	out_vcf = prex + '.vcf'
	compress(snp)
	compress(indel)
	os.system('vcf-concat %s %s > %s' %( snp_gz,indel_gz,con_vcf) )
	os.system('cat %s | vcf-sort -c > %s ' % (con_vcf,out_vcf))
	return out_vcf

@used_time
def gatk_realign_basq(dedup_bam):
	prex = dedup_bam.split('.')[0]
	tmp_name1 = prex + '.list'
	tmp_name2 = prex + '.addgrp.bam'
	tmp_name3 = prex + '.realign.bam'
	ref_dict = ref.split('.')[0] + '.dict'
	os.system('java -Xmx%s -jar %s RGPU=1 RGPL=illumina RGSM=sample1 CREATE_INDEX=True RGID=foo RGLB=ssample I=%s O=%s SORT_ORDER=coordinate'% (memory,picard_add_rg,dedup_bam,tmp_name2))
	if not os.path.exists(ref_dict):
		os.system('java -Xmx%s -jar %s R=%s O=%s' % (memory,picard_create_dict,ref,ref_dict))
	os.system('java -Xmx%s -jar %s -T RealignerTargetCreator -nt %s -R %s -I %s -o %s ' % (memory,gatk,threads,ref,tmp_name2,tmp_name1))
	os.system('java -Xmx%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s' % (memory,gatk,ref,tmp_name2,tmp_name1,tmp_name3))
	return tmp_name3
@used_time
def gatk_variant_call(ready_bam):
	prex = ready_bam.split('.')[0]
	tmp_name = prex + '.raw.vcf'
	os.system('java -Xmx%s -jar %s -T UnifiedGenotyper -nt %s -R %s -I %s -o %s   -glm BOTH ' % (memory,gatk,threads,ref,ready_bam,tmp_name))
	return tmp_name
@used_time
def filter_vcf(raw_vcf):
	prex = raw_vcf.split('.')[0]
	out = prex + '.vcf'
	os.system('java -Xmx%s -jar %s -T VariantFiltration  --variant %s -o %s -R %s --filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "snpfilter" ' % (memory,gatk,raw_vcf,out,ref))
	return out
@used_time
def move2oss(out):
	global out_dir
	out_dir = out_dir.strip('/')
	obj = '%s/%s' % (out_dir,out)
	oss.put_object_from_file(bucket,obj,out)
@used_time
def move2local(out):
	prex = out.split('.')[0]
	#os.system('mv %s %s' % (out,local_out_dir))
	#os.system('mv %s %s %s %s' % (objname1,objname2,bed,local_out_dir))
	#os.system('rm -f %s %s*' % (objname2,prex))
if __name__ == '__main__':
	start = time.time()
	print "\n\n============= basic info ===============\n\n\trun-memory:%s\n\tnum of threads :%s\n\tstart_state:%s\n\tprogram:%s\n\n\tfile1:%s\n\tfile2:%s\n\ttarget:%s" %(memory,threads,start_state,program,objname1,objname2,bed)
	if start_state == 0 and local == False and download == False and program == 'varscan' and bed == "":
		gfile1,gfile2 =  get_two_file(bucket,objname1,objname2)
		qfile1,qfile2 = uncompress_two_file(gfile1,gfile2)
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
		pileup = mpileup(dedup_bam)
		snp,indel = varscan_call(pileup)
		out = snp_indel_merge(snp,indel)
	 	move2oss(out)
	if start_state == 0 and local == False and download == False and program == 'varscan' and bed != "":
		gfile1,gfile2 =  get_two_file(bucket,objname1,objname2)
		qfile1,qfile2 = uncompress_two_file(gfile1,gfile2)
		bed =  get_one_file(bucket,bed)
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
		target_bam = intersect(dedup_bam,bed)
		pileup = mpileup(target_bam)
		snp,indel = varscan_call(pileup)
		out = snp_indel_merge(snp,indel)
		move2oss(out)
		ratio_file = compute_maped_and_target_ratio_out2file(sort_bam,target_bam)
		target_file = Compute_regions_reads(target_bam,bed)
		move2oss(ratio_file)
	 	move2oss(target_file)
	if start_state == 2 and local == False and download == False and program == 'varscan' and bed == "":
		sort_bam =  get_one_file(bucket,objname1)
		pileup = mpileup(sort_bam)
		snp,indel = varscan_call(pileup)
		out = snp_indel_merge(snp,indel)
	 	move2oss(out)
	if start_state == 2 and local == True and download == False and program == 'varscan' and bed == "":
		sort_bam = objname1
		pileup = mpileup(sort_bam)
		snp,indel = varscan_call(pileup)
		out = snp_indel_merge(snp,indel)
		move2local(out)
	if start_state == 1 and local == True and download == False and program == 'varscan' and bed == "" :
		qfile1 = objname1
		qfile2 = objname2	
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
		pileup = mpileup(dedup_bam)
		snp,indel = varscan_call(pileup)
		out = snp_indel_merge(snp,indel)
		move2local(out)
	if start_state == 0 and local == False and download == False and program == 'GATK' and bed == "":
		gfile1,gfile2 =  get_two_file(bucket,objname1,objname2)
		qfile1,qfile2 = uncompress_two_file(gfile1,gfile2)
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
	 	realign_bam = gatk_realign_basq(dedup_bam)
		raw_vcf = gatk_variant_call(realign_bam)
		filt_vcf = filter_vcf(raw_vcf)
		move2oss(filt_vcf)
	if start_state == 2 and local == False and download == False and program == 'GATK' and bed == "":
		sort_bam =  get_one_file(bucket,objname1)
	 	realign_bam = gatk_realign_basq(sort_bam)
		raw_vcf = gatk_variant_call(realign_bam)
		filt_vcf = filter_vcf(raw_vcf)
		move2oss(filt_vcf)
	if start_state == 2 and local == True and download == False and program == 'GATK' and bed == "":
		sort_bam =  objname1
	 	realign_bam = gatk_realign_basq(sort_bam)
		raw_vcf = gatk_variant_call(realign_bam)
		filt_vcf = filter_vcf(raw_vcf)
		move2local(filt_vcf)
	if start_state == 1 and local == True and download == False and program == 'GATK' and bed == "":
		qfile1 = objname1
		qfile2 = objname2	
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
	 	realign_bam = gatk_realign_basq(dedup_bam)
		raw_vcf = gatk_variant_call(realign_bam)
		filt_vcf = filter_vcf(raw_vcf)
		move2local(filt_vcf)
	if start_state == 1 and local == True and download == False and program == 'varscan' and bed != "":
		qfile1 = objname1
		qfile2 = objname2	
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
		#bed = format_bed(bed)
		target_bam = intersect(dedup_bam,bed)
		pileup = mpileup(target_bam)
		snp,indel = varscan_call(pileup)
		out = snp_indel_merge(snp,indel)
		move2local(out)
		ratio_file = compute_maped_and_target_ratio_out2file(sort_bam,target_bam)
		target_file = Compute_regions_reads(target_bam,bed)
	if start_state == 1 and local == True and download == False and program == 'GATK' and bed != "":
		qfile1 = objname1
		qfile2 = objname2	
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
	 	realign_bam = gatk_realign_basq(dedup_bam)
		target_bam = intersect(realign_bam,bed)
		raw_vcf = gatk_variant_call(target_bam)
		filt_vcf = filter_vcf(raw_vcf)
		move2local(filt_vcf)
		ratio_file = compute_maped_and_target_ratio_out2file(sort_bam,target_bam)
		target_file = Compute_regions_reads(target_bam,bed)
	if start_state == 0 and local == False and download == False and program == 'GATK' and bed != "":
		gfile1,gfile2 =  get_two_file(bucket,objname1,objname2)
		qfile1,qfile2 = uncompress_two_file(gfile1,gfile2)
		bed =  get_one_file(bucket,bed)
		sam_file = bwa_map(qfile1,qfile2)
		sort_bam = picard_sort(sam_file)
		dedup_bam= picard_dedup(sort_bam)
		target_bam = intersect(dedup_bam,bed)
	 	realign_bam = gatk_realign_basq(sort_bam)
		raw_vcf = gatk_variant_call(realign_bam)
		move2oss(raw_vcf)
		ratio_file = compute_maped_and_target_ratio_out2file(sort_bam,target_bam)
		target_file = Compute_regions_reads(target_bam,bed)
		move2oss(ratio_file)
	 	move2oss(target_file)
	end = time.time()
	used = end-start
	print "Total time is %s "  % used
