from Bio import SeqIO
import operator

fp = open("linker_lib.txt","w")
fp2 = open("rand_nucl.txt","w")
class bias_finder(object):
    records = []

    linker_len = 5
    link_dic = {}
    rand_arr = {}
    L = 5
    L_index = 51-L
    N_index = L_index - 2
    PHRED_THRESH = 2

    #Goes through a fastq file an emilinates records based on an expected
    #Phred score
    def quality_control(self):
        good_reads = (rec for rec in SeqIO.parse("example.fastq","fastq") \
                      if min(rec.letter_annotations["phred_quality"]) >= bias_finder.PHRED_THRESH)
        quality = SeqIO.write(good_reads, "good_quality.fastq","fastq")

        handle = open("good_quality.fastq","rU")
        bias_finder.records = list(SeqIO.parse(handle, "fastq"))
        handle.close()

    def link_lib(self):
        for seq_record in bias_finder.records:
            linker = str(seq_record.seq[bias_finder.L_index:51])
            rand_nucleotides = str(seq_record.seq[bias_finder.N_index:bias_finder.L_index])
            find_success = bias_finder.link_dic.get(linker, -1)
            if find_success == -1:
                bias_finder.link_dic[linker] = 1
                bias_finder.rand_arr[linker] = rand_nucleotides
            else:
                bias_finder.link_dic[linker] = bias_finder.link_dic[linker] + 1

    def print_link_lib(self):        
        sorted_x = sorted(bias_finder.link_dic.items(),key=operator.itemgetter(1),reverse=True)
        for linker, val in sorted_x:
            fp.write(linker + ": " + str(val) + "\n")
            fp2.write("Linker: " + linker + "         3' NN: " + bias_finder.rand_arr[linker] + "         Count:" + str(val) + "\n")
        fp.write("Number of unique linkers: " + str(len(sorted_x)) + "\n")
