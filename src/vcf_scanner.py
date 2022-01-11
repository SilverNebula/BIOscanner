import re
from pandas.core.series import Series
import vcf
import pandas as pd
import pyhgvs as hgvs
from pyhgvs.utils import read_transcripts
from pyfaidx import Fasta


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
        gene_detail = record.INFO['GeneDetail.refGene']
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
    vcf_reader = vcf.Reader(filename=gatk_result)
    df = pd.DataFrame(columns=['chr', 'pos', 'hgvs', 'dn', 'sg'])
    for record in vcf_reader:
        # CHROM, POS, ID, REF, ALT
        chrom = record.CHROM
        pos = record.POS
        id = record.ID
        ref = record.REF
        alt = record.ALT
        gene_detail = record.INFO['GeneDetail.refGene']
        aa_change = record.INFO['AAChange.refGene']
        # clin = record.INFO['CLNDN']
        if(record.INFO.get('CLNHGVS') is None):
            continue
        clnhgvs = record.INFO.get('CLNHGVS')
        clndn = record.INFO.get('CLNDN')
        clnsg = record.INFO.get('CLNSIG')
        if(clndn and ('Retinitis' in str(clndn))):
            df.loc[len(df)] = [chrom, pos, str(clnhgvs), clndn, clnsg]
            print('1')
            # output = '{} {}  HGVS:{} {} {}'.format(chrom,pos,clnhgvs,clndn,clnsg)
    return df


if __name__ == '__main__':
    # clinvar_result = 'D:/CODE/BIO/gene_compare/clinvar_result.txt'
    # df = pd.read_csv(clinvar_result, sep='\t')
    # df = df[['Name','Clinical significance (Last reviewed)']]
    '''
    print(df.head())
    ps = parseHGVS()
    dic = {}
    for idx,row in df.iterrows():
        name = row['Name']
        sep = str(name).find(' ')
        if sep!=-1:
            name = name[:sep]
        else:
            continue
        if name.find('(')!=-1:
            name = name[:name.find('(')] + name[name.find(')')+1:]
        # if name.find('>')==-1:
        #     continue
        if name[0]!='N':
            continue
        tran = name[:name.find(':')]
        tran = ps.gt(tran)
        if(tran == None):
            continue
        # print(name)
        # print(ps.parse(tran,name))
        dic[ps.parse(tran,name)] = row['Clinical significance (Last reviewed)']
    res= task1(dic)
    '''
    '''
    gatk_result1 = "D:/CODE/BIO/clinvar/result.vcf" #ningchunqing_FKDO210387141-1A_HMHW2DSX2_L4.hg38_multianno.vcf
    gatk_result2 = "D:/CODE/BIO/clinvar/result2.vcf" #ninghaiguang_FKDO210387142-1A_HMHW2DSX2_L4.hg38_multianno.vcf
    df1 = task2(gatk_result1)
    df2 = task2(gatk_result2)
    print(df1.head())
    df1.to_csv('csv1.csv')
    df2.to_csv('csv2.csv')
    '''
    df1 = pd.read_csv('D:/CODE/BIO/dataset/gene_compare_20211220/csv1.csv')
    df2 = pd.read_csv('D:/CODE/BIO/dataset/gene_compare_20211220/csv2.csv')
    df1 = df1[df1['hgvs'].isin(df2['hgvs'])]
    df1.to_csv('D:/CODE/BIO/dataset/gene_compare_20211220/csv_res.csv')
    # print(df2.head())
    # print(df1)
    # fout = open('clin_res.txt', 'w')
    # fout.write('\n')
    # fout.close()
    # df_res= df1[df1['hgvs']]
