#!/Users/danielence/anaconda/bin/python
#Daniel Ence
#February 4, 2019

import sys
import argparse
from cyvcf2 import VCF


def main(args):
    converter= vcf2R_qtl_CSVs(args.vcf)
    converter.set_LG_file(args.LG_file)
    converter.set_phe_file(args.phe_file)
    converter.make_CSVs(args.out)


class vcf2R_qtl_CSVs(object):

    def __init__(self,vcf):
        self.__my_vcf_file=vcf
        self.__vcf_samples = []
        self.__LG_file=""
        self.__phe_file=""
        self.__markers_not_in_map = {}

    def set_LG_file(self, LG_filename):
        self.__LG_file = LG_filename

    def set_phe_file(self, phe_filename):
        self.__phe_file = phe_filename

    def make_CSVs(self,out_prefix):
        #first write the gen_CSV
        self.__write_gen_CSV(out_prefix)
        self.__write_phe_CSV(out_prefix)

    def __write_gen_CSV(self, out_prefix):

        out_gen_filename = out_prefix + ".gen.R_qtl.csv"
        print "gen_filename is:\t" + out_gen_filename
        self.__write_R_qtl_gen_csv_header(out_gen_filename)
        self.__write_R_qtl_gen_csv_data(out_gen_filename)

    def __write_R_qtl_gen_csv_header(self,out_filename):

        out_file = open(out_filename, 'w')
        tmp_vcf_obj = VCF(self.__my_vcf_file)
        self.__vcf_samples = tmp_vcf_obj.samples

        marker_line = ""
        #make strings for all the positions
        for variant in tmp_vcf_obj:
            #if the variant is in the map, then add it to the line
            if(self.__variant_in_map(variant)):
                position = "_".join([variant.CHROM.strip().rstrip(),str(variant.start),str(variant.end)])
                marker_line = marker_line + "," + position
            #if not, then add the variant to a dictionary, indexed by scaffold name
            else:
                #print "This one was not in the map:\t" + str(variant.CHROM)
                self.__markers_not_in_map.setdefault(variant.CHROM, variant)

        print "This many markers not in the map:\t" + str(len(self.__markers_not_in_map))
        print marker_line
        out_file.write(marker_line + "\n")

        #if there is an LG_file, then get the LG's for the markers and write them out
        #if not then just print the marker line again
        if(len(self.__LG_file) > 0):
            #i know i'm actually running __make_chr_line twice to do this.
            print self.__make_LG_line()
            print self.__make_cM_line()
            out_file.write(self.__make_LG_line() + "\n")
            out_file.write(self.__make_cM_line() + "\n")
        else:
            out_file.write(marker_line)
        out_file.close()

    def __variant_in_map(self, variant):
        LG = open(self.__LG_file,'r')
        for line in LG:
            line_parts = line.strip().split("\t")
            scaffold = line_parts[13]
            if(scaffold == variant.CHROM.strip().rstrip()):
                return 1

        LG.close()
        return 0

    def __make_LG_line(self):

        LG = open(self.__LG_file, 'r')
        LG_map = {}

        for line in LG:
            line_parts = line.strip().split("\t")
            scaffold = line_parts[13]
            LG_name = line_parts[2]
            LG_pos = line_parts[3]
            if(scaffold in self.__markers_not_in_map.keys()):
                continue
            else:
                if (LG_pos == "NP"): #first trying to place NP markers in the middle of their LG
                    LG_pos = str(50) #and then will see where they shake out
                LG_map.setdefault(scaffold, {'LG_name': LG_name, 'LG_pos': LG_pos})

        tmp_VCF_obj = VCF(self.__my_vcf_file)
        chr_line = ""
        for variant in tmp_VCF_obj:
            if (variant.CHROM in LG_map):
                chr_line = chr_line + "," + LG_map[variant.CHROM.strip().rstrip()]["LG_name"]
        return chr_line

    def __make_cM_line(self):

        LG = open(self.__LG_file, 'r')
        LG_map = {}

        #comment comment

        for line in LG:
            line_parts = line.strip().split("\t")
            scaffold = line_parts[13]
            LG_name = line_parts[2]
            LG_pos = line_parts[3]
            if(scaffold in self.__markers_not_in_map.keys()):
                continue
            else:
                if (LG_pos == "NP"): #first trying to place NP markers in the middle of their LG
                    LG_pos = str(50) #and then will see where they shake out
                LG_map.setdefault(scaffold, {'LG_name': LG_name, 'LG_pos': LG_pos})

        tmp_VCF_obj = VCF(self.__my_vcf_file)
        pos_line = ""
        for variant in tmp_VCF_obj:
            if (variant.CHROM in LG_map):
                pos_line = pos_line + "," + LG_map[variant.CHROM]["LG_pos"]
        return pos_line

    def __write_R_qtl_gen_csv_data(self,out_filename):
        tmp_vcf_obj = VCF(self.__my_vcf_file)
        sample_list = tmp_vcf_obj.samples

        for i in range(len(sample_list)):
            self.__write_curr_sample_line(out_filename, sample_list[i],
                                   self.__get_sample_data(self.__my_vcf_file, i))

    def __write_curr_sample_line(self, out_filename, sample_name, sample_data):
        out_file = open(out_filename,'a')
        sample_data.insert(0,sample_name)
        out_file.write(",".join(sample_data) + "\n")

    def __get_sample_data(self, vcf_filename, sample_list_index):

        curr_sample_data = []
        for variant in VCF(vcf_filename):
            if(variant.CHROM not in self.__markers_not_in_map.keys()):
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

    def __write_phe_CSV(self, out_prefix):
        out_phe_filename = out_prefix + ".phe.R_qtl.csv"
        out_phe_file = open(out_phe_filename, 'w')

        #write out the data from the phe file and from the
        #markers_to_map in the phe csv file

        phe_file_lines = []
        if(self.__phe_file):
            phe_file_lines = self.__get_phenotype_data(out_phe_file)
        else:
            phe_file_lines.append(["ID"])
            for sample in self.__vcf_samples:
                phe_file_lines.append([sample])

        print phe_file_lines
        phe_file_lines = self.__write_markers_not_in_map(phe_file_lines)
        print phe_file_lines

        for line in phe_file_lines:
            print ",".join(line)
            out_phe_file.write(",".join(line) + "\n")
        out_phe_file.close()

    def __get_phenotype_data(self, out_file):

        tmp_vcf_sample_map = {}
        for sample in self.__vcf_samples:
            tmp_vcf_sample_map.setdefault(sample,1)

        phe_file_lines = []
        phe_file = open(self.__phe_file)
        for line in phe_file:
            line_parts = line.split("\t")

            #print str((line_parts[0] == "SAMPLE" or line_parts[0] in tmp_vcf_sample_map))
            if(line_parts[0] == "ID" or line_parts[0] in tmp_vcf_sample_map):
                phe_file_lines.append(line.strip().split("\t"))
        return phe_file_lines

    def __write_markers_not_in_map(self, phe_file_lines):

        for variant in self.__markers_not_in_map:
            variant_obj = self.__markers_not_in_map[variant]

            sample_counter = 0
            for GT in variant_obj.gt_types:
                if(sample_counter == 0):
                    #code the GT by the number of alternate alleles
                    #append the scaffold names (as the phenotype name)
                    #to the first line in the "phe_file_lines" structure
                    phe_file_lines[sample_counter].append(variant_obj.CHROM)
                    sample_counter = sample_counter + 1

                if(GT == 2):
                    GT == "NA"
                elif(GT == 3):
                    GT == 2
                phe_file_lines[sample_counter].append(str(GT))
                sample_counter = sample_counter + 1

        return  phe_file_lines

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",type=str,help="vcf file to convert for R/qtl")
    parser.add_argument("--LG_file",type=str,help="two-column tab-delim file with \"scaffold\\tLG\"")
    parser.add_argument("--phe_file", type=str,help="input file with sample phenotypes")
    parser.add_argument("--out",type=str,help="output prefix for the gen and phe files. Files will end with \".R_qtl.csv\"")
    args = parser.parse_args()
    main(args)
