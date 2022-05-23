import re
import numpy
from pandas.core.series import Series
import vcf
import pandas as pd
import pyhgvs as hgvs
from pyhgvs.utils import read_transcripts
from pyfaidx import Fasta
import argparse


class parseHGVS():
    def __init__(self):
        # Read genome sequence using pyfaidx.
        self.genome = Fasta('D:/CODE/BIO/hg38.fa')

        # Read RefSeq transcripts into a python dict.
        with open('./genes.refGene') as infile:
            self.transcripts = read_transcripts(infile)

        # Provide a callback for fetching a transcript by its name.
        def get_transcript(name):
            return self.transcripts.get(name)
        self.gt = get_transcript
        pass

    def parse(self, tran, name):
        # Parse the HGVS name into genomic coordinates and alleles.
        chrom, offset, ref, alt = hgvs.parse_hgvs_name(
            name, self.genome, transcript=tran, get_transcript=self.gt)
        # Returns variant in VCF style: ('chr11', 17496508, 'T', 'C')
        return chrom, offset, ref, alt


def task1(dic: dict):
    gatk_result = "D:/CODE/BIO/gene_compare/ninghaiguang_FKDO210387142-1A_HMHW2DSX2_L4.hg38_multianno.vcf"
    vcf_reader = vcf.Reader(filename=gatk_result)
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        ref = record.REF
        alt = record.ALT
        gene_detail = record.INFO['GeneDetail.GENEINFO']
        aa_change = record.INFO['AAChange.refGene']
        res = (chrom, pos, ref, str(alt[0]))
        # print(res)
        res = dic.get(res)
        if(res):
            print(res)

    '''
    fout = open('result.txt', 'w')
    for name in st:
        fout.write(name)
        fout.write('\n')
    fout.close()
    '''


def task2(gatk_result):
    vcf_reader = vcf.Reader(filename=gatk_result,encoding='utf8')
    df = pd.DataFrame(columns=['chr', 'pos','geneinfo', 'hgvs', 'dn', 'sg'])
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        ref = record.REF
        alt = record.ALT
        gene_detail = record.INFO.get('GENEINFO')
        # aa_change = record.INFO['AAChange.refGene']
        # clin = record.INFO['CLNDN']
        if(record.INFO.get('CLNHGVS') is None):
            continue
        clnhgvs = record.INFO.get('CLNHGVS')
        clndn = record.INFO.get('CLNDN')
        clnsg = record.INFO.get('CLNSIG')
        if(('MT-ND6' in str(gene_detail))):
            df.loc[len(df)] = [chrom, pos,gene_detail, str(clnhgvs), clndn, clnsg]
            print('1')
            # output = '{} {}  HGVS:{} {} {}'.format(chrom,pos,clnhgvs,clndn,clnsg)
    return df


def find_splice_variant(vcf_db):
    vcf_reader = vcf.Reader(filename=vcf_db,encoding='utf8')
    fout = open('result_false.vcf','w')
    vcf_writer = vcf.Writer(fout,vcf_reader)
    output =[]
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        ref = record.REF
        if len(ref)!=1:
            continue
        alt = record.ALT
        if (alt[0] is None) or (len(alt[0]))!=1:
            continue
        gene_detail = record.INFO.get('GENEINFO')
        clnhgvs = record.INFO.get('CLNHGVS')
        clndn = record.INFO.get('CLNDN')
        clnsg = record.INFO.get('CLNSIG')
        if clnsg and (clnsg[0] not in 'Pathogenic'):
            continue
        vtype = record.INFO.get('MC')
        if vtype is None:
            continue
        record.INFO = {}
        record.add_info('CLNSIG',clnsg)
        record.add_info('CLNDN',clndn)
        record.add_info('MC',vtype)
        record.add_info('CLNHGVS',clnhgvs)
        ### change filter here
        if (('splice' not in vtype[0]) and (('missense_variant' in vtype[0])or('synonymous_variant' in vtype[0]))):
            output.append(record)
        ###
    idx = numpy.random.permutation(len(output)).tolist()
    for i in range(8000):
        vcf_writer.write_record(output[idx[i]])
    vcf_writer.close()
    

def match_splice_variant(vcf_db):
    vcf_reader = vcf.Reader(filename=vcf_db,encoding='utf8')
    fout = open('result_false_filtered','w')
    hgvs_dic = {}
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        clnhgvs = record.INFO.get('CLNHGVS')
        hgvs_dic[(chrom,pos,id)] = clnhgvs
    vcf_reader = vcf.Reader(filename='result_false.vcf',encoding='utf8')
    vcf_writer = vcf.Writer(fout,vcf_reader)
    output =[]
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        ref = record.REF
        hgvs = hgvs_dic.get((chrom,pos,id))
        if hgvs is None:
            continue
        if len(ref)!=1:
            continue
        alt = record.ALT
        if (alt[0] is None) or (len(alt[0]))!=1:
            continue
        gene_detail = record.INFO.get('GENEINFO')
        clndn = record.INFO.get('CLNDN')
        clnsg = record.INFO.get('CLNSIG')
        vtype = record.INFO.get('MC')
        record.INFO = {}
        record.add_info('CLNSIG',clnsg)
        record.add_info('CLNDN',clndn)
        record.add_info('MC',vtype)
        record.add_info('CLNHGVS',hgvs)
        ### change filter here
        output.append(record)
        ###
    idx = numpy.random.permutation(len(output)).tolist()
    for i in range(len(output)):
        vcf_writer.write_record(output[idx[i]])
    vcf_writer.close()
    pass

def calc_predication(result_vcf):
    vcf_reader = vcf.Reader(filename=result_vcf,encoding='utf8')
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        ref = record.REF
        alt = record.ALT
        score = ''
        for k in record.INFO.keys():
            if record.INFO[k] == True:
                score = k
                break
                # print(k)
        score = score.split(sep='|')[2:6]
        score = [float(x) for x in score]
        print(score)

if __name__ == '__main__':
    '''
    parser = argparse.ArgumentParser(description="training setting")
    parser.add_argument('targetfile',type=str,help='target .vcf file')
    args = parser.parse_args()
    gatk_result2 = args.targetfile
    
    df2 = task2(gatk_result2)
    print(df2.head())
    df2.to_csv('csv2.csv')
    '''
    path = r'D:\CODE\BIO\dataset\clinvar_dataset_2022\grch37\clinvar_20220507.vcf'
    match_splice_variant(path)
    # path2 = 'splice_result.vcf'
    # calc_predication(path2)
