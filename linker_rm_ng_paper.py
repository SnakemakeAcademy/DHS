# -*- coding: utf-8 -*-

import argparse
import gzip
import datetime
# Import pairwise2 module
#from Bio import pairwise2


def remove_linker(options):

    total = 0
    linker_seq = 0
    linker_not_correct_loc = 0
    linker_non = 0
    output1 = options.prefix + "_left.fq"
    output2 = options.prefix + "_right.fq"
    output_log = options.prefix + ".log"
    log = open(output_log, 'w')
    log.write("Total\tLinker_hit\tLinker_hit_notgood\tLinker_non\n")
    flag = 0
    with open(options.read) as f1, open(output1, 'w') as output_f1, open(output2, 'w') as output_f2:
        for id1 in f1:
            flag = flag + 1

            if flag % 100000 == 0:
                currentDT = datetime.datetime.now()
                #log.write(str(currentDT) + ": " + str(flag) + "\n")
                #log.flush()

            seq1 = next(f1)
            plus1 = next(f1)
            qual1 = next(f1)
            #id1 = id1.decode('utf-8')
            #seq1 = seq1.decode('utf-8')
            seq1 = seq1.rstrip()
            qual1 = qual1.rstrip()
            if len(seq1) < 75 or len(seq1) > 81:
                continue
            #alignments = pairwise2.align.localms(seq1, options.linker, 2, -1, -0.5, -0.1)
            index = seq1.find(options.linker)
            total = total + 1
            if index > -1:
                if(index >=17 and  index <= 21):
                    output_f1.write(id1)
                    seq_r1 = seq1[0:index]
                    output_f1.write(seq_r1 + "\n")
                    output_f1.write("+\n")
                    qual1_r1 = qual1[0:index]
                    output_f1.write(qual1_r1 + "\n")
                    output_f2.write(id1)
                    linker_end = index + len(options.linker)
                    seq_r2 = seq1[linker_end:]
                    output_f2.write(seq_r2 + "\n")
                    output_f2.write("+\n")
                    qual1_r2 = qual1[linker_end:]
                    output_f2.write(qual1_r2 + "\n")
                    linker_seq = linker_seq + 1
                else:
                    linker_not_correct_loc = linker_not_correct_loc + 1
            else:
                linker_non = linker_non + 1
    log.write("{total_num}\t{ls}\t{l_not_enough}\t{l_non}\n".format(total_num = total, ls = linker_seq, l_not_enough = linker_not_correct_loc, l_non = linker_non))
    log.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--read", help = "Read",required=True)
    parser.add_argument("-l", "--linker", help = "Linker sequences", default = "TAGTCGGAGAACCAGTAGCTAGCTACTGGTTCTCCGAC")
    parser.add_argument("-p", "--prefix", help = "Prefix", default="linker_rm")
    parser.add_argument("-c", "--cutoff", help = "Cutoff", default=75, type=int)
    options = parser.parse_args()
    remove_linker(options)
