import numpy as np
import esm
import io
import argparse


class Seq:
    def __init__(self, sequence, pred):
        self.sequence = sequence
        self.base = np.array(self.sequence, dtype=np.str)
        self.value = np.array(pred)
        self.len = len(sequence)
        ############################
        # Configuration
        self.open_threshold = 0.6
        self.open_len_threshold = 8
        # ratio of 'open' nucleotide should between [lower_bound,upper_bound]
        self.upper_bound = 1.0
        self.lower_bound = 0.3
        # window size should between [min,max)
        self.max_window = 21
        self.min_window = 20
        # If there are windows with start(or end) index difference smaller than resolution, only one of them will be count
        # e.g : If resolution==3 and there are 4 windows (1,3),(1,4),(1,5),(1,6), only (1,3) and (1,6) will be outputed
        self.resolution = 1
        # value of Tm temperature should be equal or above threshold
        self.Tm_threshold = 48

        # value of (G+C)/windowsize should between [lowerbound,upperbound]
        self.CG_windowsize = 20
        self.CG_upperbound = 0.6
        self.CG_lowerbound = 0.4
        #
        self.forbidden_seq = ['GGGG']
        self.positive_seq = ['CCAC', 'TCCC', 'ACTC', 'GCCA', 'CTCT']
        self.negative_seq = ['GGGG', 'ACTG', 'TAA', 'CCGG', 'AAA']
        ############################
        self.init_automaton()
        pass

    def init_automaton(self):
        self.forbid_index = esm.Index()
        for word in self.forbidden_seq:
            self.forbid_index.enter(word)
        self.forbid_index.fix()

        self.positive_index = esm.Index()
        for word in self.positive_seq:
            self.positive_index.enter(word)
        self.positive_index.fix()

        self.negative_index = esm.Index()
        for word in self.negative_seq:
            self.negative_index.enter(word)
        self.negative_index.fix()
        pass

    def get_sub_str(self, st, ed):
        return ''.join(reversed(self.sequence[st:ed]))

    def get_len(self):
        return self.len

    def calc_contiguous_opening_len(self, st, ed):
        len = 0
        res = 0
        for i in range(st, ed):
            if(self.value[i] >= self.open_threshold):
                len += 1
            else:
                len = 0
            res = max(res, len)
        return res


def calcGC(seq: Seq, st, ed):
    # examine G+C ratio and Tm
    threshold2 = seq.CG_windowsize * seq.CG_upperbound
    threshold1 = seq.CG_windowsize * seq.CG_lowerbound
    cnt = 0
    for i in range(st, ed):
        # count the element that moves into window
        cnt += (seq.base[i] == 'C') + (seq.base[i] == 'G')
        # count the element that just move out of window
        if(i-seq.CG_windowsize >= st):
            cnt -= (seq.base[i-seq.CG_windowsize] == 'C') + \
                (seq.base[i-seq.CG_windowsize] == 'G')
        #
        if(i-seq.CG_windowsize+1 >= st):  # a window
            if((cnt < threshold1) or (cnt > threshold2)):
                return False
            Tm = 4 * cnt + 2 * (seq.CG_windowsize - cnt)
            if(Tm < seq.Tm_threshold):
                return False
    return True


def calc_special_subseq(seq: Seq, st, ed):
    query_str = ''.join(seq.sequence[st:ed])
    forbidden = seq.negative_index.query(query_str)
    positive = seq.positive_index.query(query_str)
    negative = seq.negative_index.query(query_str)
    return len(forbidden), len(positive), len(negative)


def calc(seq: Seq, output_path):
    fp = None
    if output_path:
        fp = open(output_path, 'w')
    ans_set = set()
    for i in range(seq.len):
        for j in range(i + 1, min(i + seq.max_window, seq.len)):
            open_cnt = np.sum(seq.value[i:j] >= seq.open_threshold)
            close_cnt = np.sum(seq.value[i:j] < seq.open_threshold)
            ratio = open_cnt / (j-i)
            if((ratio < seq.lower_bound) or (ratio > seq.upper_bound)):
                continue
            contiguous_len = seq.calc_contiguous_opening_len(i, j)
            if(contiguous_len < seq.open_len_threshold):
                continue
            if((j-i) < seq.min_window):
                continue
            # if(calcGC(seq, i, j) == False):
            #     continue
            fb, posi, nega = calc_special_subseq(seq, i, j)
            if(fb > 0):
                continue
            if((i//seq.resolution, j//seq.resolution) not in ans_set):
                ans_set.add((i//seq.resolution, j//seq.resolution))
                '''
                output_str = '{} - {} \t len: {} \t positive count:{} \t negative count:{} \t content:{}'.format(
                    i+1, j+1, (j-i), posi, nega, seq.get_sub_str(i, j))
                print(output_str)
                '''
                # i+1 because array index starts from 0
                fp.write('{}-{},'.format(i+1, j+1))
                # fp.write('\n')
        pass
    print("Total number:")
    print(len(ans_set))
    fp.write("\nTotal number: {}\n".format(len(ans_set)))
    pass


if __name__ == '__main__':
    char_map = {'A': 'T',
                'G': 'C',
                'C': 'G',
                'U': 'A'}

    print('start')
    seq = None
    #
    sscount_file = './sscount.txt'
    output_file = './out.txt'
    #
    with open(sscount_file, 'r') as fp:
        cnt = fp.readline()
        data = fp.readlines()
        print(int(cnt))
        value = []
        sequence = []
        for line in data:
            sp = line.split()
            value.append(int(sp[1])/int(cnt))
            sequence.append(char_map[sp[2]])
        seq = Seq(sequence, value)
    calc(seq, output_file)
