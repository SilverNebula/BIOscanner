import xml.etree.cElementTree as ET
import os


def parse_name(name):
    res = str(name)
    pos = str.find(name, ':')
    if pos > -1:
        res = name[pos+1:] + ':' + name[:pos]
    return res


def get_association_info(association_node):
    gene_node = association_node.find('Gene')
    name = gene_node.find('Symbol').text
    return name


def parse_association_list(list_node):
    genes = []
    for node in list_node.iter('DisorderGeneAssociation'):
        genes.append(get_association_info(node))
    return genes


class Disorder:
    def __init__(self, node=None):
        if node:
            name = node.find('Name').text
            orpha_code = node.find('OrphaCode').text
            expert_link = node.find('ExpertLink').text
            genes_array = parse_association_list(
                node.find('DisorderGeneAssociationList'))
            #
            self.name = name
            self.orpha_code = int(orpha_code)
            self.expert_link = expert_link
            self.genes_array = genes_array
            self.disorder_type = node.find('DisorderType').find('Name').text
            self.disorder_group = node.find('DisorderGroup').find('Name').text
            pass
        pass

    def __lt__(self, other):
        return int(self.orpha_code) < int(other.orpha_code)
    pass


def get_disorders(root):
    res = []
    for disorder in root:
        temp = Disorder(disorder)
        res.append(temp)
    res.sort()
    return res


def generate_disorder_txt(fout, disorder_list):
    # output
    for x in head:
        fout.write('{}: {}\n'.format(x, head[x]))
    fout.write('-----------------------\n')
    for disorder in disorder_list:
        genes_array_str = ",".join(str(x) for x in disorder.genes_array)
        fout.write('Disorder Name: {}\n - orpha_code: {}\n - expert_link: {}\n'.format(
            parse_name(disorder.name), disorder.orpha_code, disorder.expert_link))
        fout.write(' - disorder_type: {}\n - disorder_group: {}\n'.format(
            disorder.disorder_type, disorder.disorder_group))
        fout.write(' - associated genes: {}\n'.format(genes_array_str))
        fout.write('\n')


def search_gene_by_disorder(fout, disorder_list, name):
    for disorder in disorder_list:
        if name.lower() in disorder.name.lower():
            genes_array_str = "\n".join(str(x) for x in disorder.genes_array)
            fout.write(genes_array_str)
            fout.write('\n')
    return


def search_disorder_by_gene(fout, disorder_list, gene_file):
    gene_file = "D:/CODE/BIO/gene_compare/result.txt"
    fin = open(gene_file, 'r')
    gene_set = set()
    #
    line = fin.readline()
    while(line):
        name = line.strip()
        gene_set.add(name)
        line = fin.readline()


if __name__ == '__main__':
    fout = open('data.txt', "w", encoding='UTF-8')
    # tree = ET.parse('./en_product6.xml')
    tree = ET.parse('./original.xml')
    root = tree.getroot()
    # process head
    head = root.attrib
    print(head)
    # process data
    root = root.find('DisorderList')
    disorder_list = get_disorders(root)
    # generate_disorder_txt(fout,disorder_list)
    search_gene_by_disorder(fout, disorder_list, 'Retinitis Pigmentosa')
    fout.close()
    pass
