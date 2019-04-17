#!/Users/danielence/anaconda/bin/python
#Daniel Ence
#March 1, 2019

import sys
import argparse
from cyvcf2 import VCF


def main(args):

    converter = vcf2onemap(args.vcf)
    converter.make_onemap_raw(args.out)

class vcf2onemap(object):

    def __init__(self,vcf):
        self.__my_vcf_file = vcf

    def make_onemap_raw(self, out_prefix):
        outfile = open(out_prefix + ".onemap.raw",'w')
        self.__write_onemap_header(outfile)
        self.__write_onemap_data(outfile)

    def __write_onemap_header(self,out_handle):

        vcf_file = VCF(self.__my_vcf_file)
        variant_counter = 0
        sample_counter = 0
        for variant in vcf_file:
            if(sample_counter == 0):
                sample_counter = len(variant.gt_types)
            variant_counter = variant_counter + 1

        header_line = str(sample_counter) + " " + str(variant_counter)
        out_handle.write(header_line + "\n")

    def __write_onemap_data(self, out_handle):

        vcf_file = VCF(self.__my_vcf_file)

        for variant in vcf_file:
            out_handle.write(self.__make_onemapraw_line(variant) + "\n")

    def __make_onemapraw_line(self, variant_obj):

        onemapraw_line = "*"
        onemapraw_line = onemapraw_line + self.__make_marker_name(variant_obj)

        marker_type = self.__determine_marker_type(variant_obj)

        onemapraw_line = onemapraw_line + " " + marker_type

        first_marker = True
        for gt_type in variant_obj.gt_types:

            allele = ""
            if(gt_type == 3):
                allele = "b"
            elif(gt_type == 2):
                allele = "-"
            elif(gt_type == 1):
                allele = "ab"
            elif(gt_type == 0):
                allele = "a"

            if first_marker == True:
                onemapraw_line = onemapraw_line + " " + allele
                first_marker = False
            else:
                onemapraw_line = onemapraw_line + "," + allele

        return onemapraw_line

    def __make_marker_name(self,variant_obj):
        return "_".join((variant_obj.CHROM,str(variant_obj.start)))

    def __determine_marker_type(self, variant_obj):

        HET_FOUND = False

        for GT in variant_obj.gt_types:
            if(GT == 1):
                HET_FOUND = True

        if HET_FOUND == True:
            return "B3.7"
        else:
            return "D1.11"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, help="vcf file to convert for onemap raw format")
    parser.add_argument("--out", type=str,help="output prefix for the \"onemap.raw\" file that will be made")
    args = parser.parse_args()
    main(args)