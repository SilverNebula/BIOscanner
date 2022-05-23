import requests
from lxml import etree
import numpy as np
import vcf
data = {
    'organism':'human',
    'which':'acceptors',
    'reverse':'no',
    'min_donor':0.4,
    'min_acc':0.4,
    'text':'AAAAAAAAAAAAAAAAAAAAAAAG'
}
proxies = {
    'http': '127.0.0.1:10809/',
    'https': '127.0.0.1:10809/'
}
cookie ={
    '_ga':'GA1.2.1929991677.1652851777',
    '_gid':'GA1.2.1922660584.1652943850',
}
header = {
    'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.67 Safari/537.36'
}


target_vcf = 'result_false_filtered.vcf'
vcf_reader = vcf.Reader(filename=target_vcf,encoding='utf8')
fout = open('result_false_final.vcf','w')
vcf_writer = vcf.Writer(fout,vcf_reader)
for record in vcf_reader:
    hgvs = record.INFO.get('CLNHGVS')
    id = record.ID
    url = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + str(id)
    print(hgvs)
    print(id)
    print(url)
    resp = requests.get(url=url,proxies=proxies,verify=False,timeout=5)
    tree = etree.HTML(resp.text)
    nodes = tree.xpath('//div[@class="full_comment collapsed"]/text()')
    splice_flag = 0
    for node in nodes:
        if node and ('splic' in node):
            splice_flag = True
            print(node)
            break
    if splice_flag==0:
        vcf_writer.write_record(record)
vcf_writer.close()