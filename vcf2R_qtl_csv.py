#!/Users/danielence/anaconda/bin/python
#Daniel Ence
#February 4, 2019

import sys
import argparse
from cyvcf2 import VCF


def main(args):
    converter= vcf2R_qtl(args.vcf)
    converter.set_LG_file(args.LG_file)
    converter.make_csv(args.out)


class vcf2R_qtl(object):

    def __init__(self,vcf):
        self.__my_vcf_file=vcf
        self.__LG_file=""

    def set_LG_file(self, LG_filename):
        self.__LG_file = LG_filename
        print "LG_file is:\t"
        print self.__LG_file

    def make_csv(self,out_prefix):

        out_filename = out_prefix + ".R_qtl.csv"
        out_file =  open(out_filename,'w')

        self.write_R_qtl_csv_header(out_file)
        self.write_R_qtl_csv_data(out_file)

    def write_R_qtl_csv_data(self,out_file):
        tmp_vcf_obj = VCF(self.__my_vcf_file)
        sample_list = tmp_vcf_obj.samples

        for i in range(len(sample_list)):
            self.__write_curr_sample_line(out_file, sample_list[i],
                                   self.__get_sample_data(self.__my_vcf_file, i))

    def __write_curr_sample_line(self, out_file, sample_name, sample_data):
        sample_data.insert(0,sample_name)
        out_file.write(",".join(sample_data) + "\n")

    def __get_sample_data(self, vcf_filename, sample_list_index):

        curr_sample_data = []
        for variant in VCF(vcf_filename):
            curr_sample_gt = variant.gt_types[sample_list_index]
            if(curr_sample_gt == 0):
                curr_sample_data.append("A")
            elif(curr_sample_gt == 1):
                curr_sample_data.append("AB")
            elif(curr_sample_gt == 2):
                curr_sample_data.append("NA")
            elif(curr_sample_gt == 3):
                curr_sample_data.append("B")
        return curr_sample_data

    def write_R_qtl_csv_header(self,out_file):

        tmp_vcf_obj = VCF(self.__my_vcf_file)
        marker_line = ""
        #make strings for all the positions
        for variant in tmp_vcf_obj:
            position = "_".join([variant.CHROM,str(variant.start),str(variant.end)])
            marker_line = marker_line + "," + position

        #if there is an LG_file, then get the LG's for the markers and writ them out
        #if not then just print the marker line again
        out_file.write(marker_line + "\n")
        if(len(self.__LG_file) > 0):
            #i know i'm actually running __make_chr_line twice to do this.
            out_file.write(self.__make_chr_line()[0] + "\n")
            out_file.write(self.__make_chr_line()[1] + "\n")
        else:
            out_file.write(marker_line)

    def __make_chr_line(self):
        LG = open(self.__LG_file,'r')
        LG_map = {}
        for line in LG:
            line_parts = line.strip().split("\t")
            scaffold = line_parts[0]
            LG_name = line_parts[1]
            SNP_type = line_parts[2]
            LG_pos = line_parts[3]
            if(LG_pos == "NP"):
                LG_pos = str(50)
            LG_map.setdefault(scaffold,{'LG_name': LG_name,'LG_pos': LG_pos})

        tmp_VCF_obj = VCF(self.__my_vcf_file)
        chr_line = ""
        pos_line = ""
        for variant in tmp_VCF_obj:

            if(variant.CHROM in LG_map):
                chr_line = chr_line + "," + LG_map[variant.CHROM.strip().rstrip()]["LG_name"]
                pos_line = pos_line + "," + LG_map[variant.CHROM.strip().rstrip()]["LG_pos"]
            else:
                chr_line = chr_line + "," + "no_LG"
                pos_line = pos_line + "," + "50"

        return (chr_line,pos_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",type=str,help="vcf file to convert for R/qtl")
    parser.add_argument("--LG_file",type=str,help="two-column tab-delim file with \"scaffold\\tLG\"")
    parser.add_argument("--out",type=str,help="output prefix for the gen and phe files. Files will end with \".R_qtl.csv\"")
    args = parser.parse_args()
    main(args)
