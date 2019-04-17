#!/Users/danielence/anaconda/bin/python
#Daniel Ence
#May 31, 2018

import argparse
from cyvcf2 import VCF

def main(args):

    my_filemaker = locus_genotype_filemaker(args.vcf,args.name,args.popt,args.sample_list,args.interval_file, args.phenotype_file)
    if(my_filemaker.get_filename() == ""):
        print "the filename is empty"
    else:
        print "the filename is:\t" + my_filemaker.get_filename()
    print "poptype is:\t" + my_filemaker.get_poptype()
    print "name is:\t" + my_filemaker.get_name()

    my_filemaker.set_LG_group_file(args.LG_map_file)
    my_filemaker.set_subset_file(args.subset_samples)
    my_filemaker.make_joinmap_locus_genotype_file()

class locus_genotype_filemaker(object):

    #always use an empty constructor to just start out
    def __init__(self):
        self.__my_vcf_file=""
        self.__my_name=""
        self.__my_poptype=""
        self.__my_target_indvs=[]
        self.__my_target_intervals=[]
        self.__LG_group_file = ""
        self.__subset_file = ""
        self.__phenotype_map = ""

    def __init__(self,vcf_filename,name,poptype,group_arg,intervals_file,phenotype_file):
        self.__my_vcf_file=vcf_filename
        self.__my_name=""
        self.__LG_group_file = ""
        self.__subset_file = ""

        if(self.__is_valid_name(name)):
            self.__my_name = name
        self.__my_poptype=""
        if(self.__is_valid_poptype(poptype) == 1):
            self.__my_poptype=poptype

        self.__my_target_indvs = []
        if(group_arg != None):
            self.__my_target_indvs = self.setup_group_member_list(group_arg)

        self.__my_target_intervals = []
        if(intervals_file != None):
            self.__my_target_intervals = self.__setup_target_intervals(intervals_file)

        self.__phenotype_map = {}
        if(phenotype_file != None):
            self.__phenotype_map = self.__setup_phenotype_map(phenotype_file)

    def __setup_phenotype_map(self,phenotype_file):
        #expect a two column, tab-delim file.
        #column 1 is the ID (or longer ID)
        #column 2 is the pheontype (binary)
        pheno_file_handle = open(phenotype_file,'r')
        pheno_map = {}
        for line in pheno_file_handle:
            line_parts=line.strip().split("\t")
            pheno_map[line_parts[0]] = line_parts[1]
        pheno_file_handle.close()

        #print pheno_map
        return pheno_map

    def get_name(self):
        return self.__my_name

    def get_filename(self):
        return self.__my_vcf_file

    def get_poptype(self):
        return self.__my_poptype

    def __is_valid_name(self,poss_name):

        if(poss_name.replace(' ',"_") == poss_name):
            return 1
        else:
            return 0

    def __is_valid_poptype(self,poss_poptype):
        #I think it is ok to restrict the scope of this dictionary to just this method
        valid_poptypes_dict={"BC1":1,"F2":1,"RIx":1,"DH":1,"DH1":1,"DH2":1,"HAP":1,"HAP1":1,"CP":1,"BCpxFy":1,"IMxFy":1}

        if(valid_poptypes_dict.get(poss_poptype) ==1):
            return 1
        else:
            return 0

#Need to put this method in a utils for all the vcf converter
    def setup_group_member_list(self, group_arg_string):

        group_list = self.get_list_from_file(group_arg_string)

        group_members_hash = {}
        for tmp_member in group_list:
            group_members_hash.setdefault(tmp_member, 1)

        group_members_list = []

        vcf = VCF(self.__my_vcf_file)

        for sample in vcf.samples:
            if (group_members_hash.get(sample) == 1):
                group_members_list.append(1)
            else:
                group_members_list.append(0)

        return group_members_list

    def get_list_from_file(self, filename):

        members_list = []
        for line in open(filename, 'r'):
            line = line.strip()
            # allow for multiple columns in the file, but assume the first one is the sample name
            parts = line.split("\t")
            members_list.append(parts[0])
        return members_list

    #this filter hasn't been implemented yet, DE 05/31/2018
    def __setup_target_intervals(self, intervals_file):

        #expecting a properly formatted bedfile. We'll only use the first three fields
        intervals = open(intervals_file,'r')
        for line in intervals:
            parts = line.strip().split("\t")
            self.__my_target_intervals.append((parts[0],parts[1],parts[2]))

    def get_indvs_count(self):
        if(len(self.__my_target_indvs) > 0):
            counter = 1
            counter_zero = 1
            for i in range(len(self.__my_target_indvs)):
                    if(self.__my_target_indvs[i] == 1):
                            counter = counter + 1
                    else:
                            counter_zero = counter_zero + 1

            return counter
        else:
            tmp_vcf = VCF(self.__my_vcf_file)
            return len(tmp_vcf.samples)

    def get_intervals_count(self):
        if(len(self.__phenotype_map) == 0):
            if(len(self.__my_target_intervals) > 0):
                return self.__my_target_intervals.len()
            else:
                tmp_vcf = VCF(self.__my_vcf_file)
                counter = 1
                for variant in tmp_vcf:
                    counter = counter + 1
                return counter
        else:
            if(len(self.__my_target_intervals) > 0):
                return self.__my_target_intervals.len() + 1
            else:
                tmp_vcf = VCF(self.__my_vcf_file)
                counter = 2
                for variant in tmp_vcf:
                    counter = counter + 1
                return counter

    def make_joinmap_locus_genotype_file(self):

        locus_genotype_file = self.__my_name + ".jm4.with_pheno.txt"
        print "genotype file is:\t" + locus_genotype_file
        outfile = open(locus_genotype_file,'w')

        self.__write_locus_genotype_file_header(outfile)
        self.__write_locus_genotype_file_data(outfile)
        if(len(self.__phenotype_map.keys()) > 0):
            self.__write_phenotype_data(outfile)

        self.__write_sample_names(outfile)

        outfile.close()

    def __write_sample_names(self, outfile):

        outfile.write("\nindividual names:\n")
        tmp_vcf = VCF(self.__my_vcf_file)

        if(len(self.__my_target_indvs) > 0):
            for i in range(len(self.__my_target_indvs)):
                if(self.__my_target_indvs[i] == 1):
                    outfile.write(tmp_vcf.samples[i] + "\n")
        else:
            tmp_vcf = VCF(self.__my_vcf_file)
            for sample in tmp_vcf.samples:
                outfile.write(sample+ "\n")

    def __write_locus_genotype_file_header(self,outfile_handle):
        outfile_handle.write("name = " + self.__my_name + "\n")
        outfile_handle.write("popt = " + self.__my_poptype + "\n")

        if (len(self.__my_target_intervals) > 0):
            outfile_handle.write("nloc = " + str(len(self.__my_target_intervals)) + "\n")
        else:
            outfile_handle.write("nloc = " + str(self.get_intervals_count()) + "\n")

        if (len(self.__my_target_indvs) > 0):
            outfile_handle.write("nind = " + str(self.get_indvs_count()) + "\n")
        else:
            outfile_handle.write("nind = " + str(self.get_indvs_count()) + "\n")

    def __write_locus_genotype_file_data(self, outfile):

        tmp_vcf = VCF(self.__my_vcf_file)
        counter = 0
        variant_counter_filehandle = open("variant_counter_map.txt", 'w')
        for variant in tmp_vcf:
            outfile.write(self.__make_locus_genotype_file_data_line(variant, counter, variant_counter_filehandle) + "\n")
            counter = counter + 1
        variant_counter_filehandle.close()

    def __write_phenotype_data(self,outfile):
        #assume binary phenotypes for now
        #pheno_1=""
        #pheno_2=""
        #for val in list(self.__phenotype_map.values()):
        #    if(len(pheno_1) == 0):
        #        pheno_1 = val
        #    else:
        #        if(val != pheno_1):
        #            pheno_2 = val

        #need to match the order of the samples in the vcf header
        tmp_vcf = VCF(self.__my_vcf_file)
        tmp_pheno_list = []
        for sample in tmp_vcf.samples:
            if(self.__phenotype_map.get(sample) != None):
                if(self.__phenotype_map.get(sample) == "gall"):
                    tmp_pheno_list.append("ll")
                else:
                    tmp_pheno_list.append("lm")
        fields_list = ["phenotype", "<lmxll>","{1-}","(ll,lm)"] + tmp_pheno_list
        line = " ".join(fields_list)
        outfile.write(line + "\n")

    def __make_locus_genotype_file_data_line(self,variant, counter, counter_filehandle):
            locus_name = self.__make_locus_name("SNP",counter,variant)
            counter_filehandle_line = "\t".join((locus_name,variant.CHROM,str(variant.start),str(variant.end)))
            #print "locus_name is:\t" + locus_name
            #print variant.CHROM
            #print variant.start
            #print variant.end
            counter_filehandle.write(counter_filehandle_line + "\n")

            #hardcoding this in for my project even though its bad practice
            #could make this a parameter for the script maybe
            segregation,clas,phase = self.__determine_segregation(variant.gt_types)

            fields_list = []
            if(segregation != None):
                fields_list = [locus_name,segregation,phase,self.__clas_to_string(clas)]
            else:
                fields_list = [locus_name,phase,self.__clas_to_string(clas)]

            for i in range(len(variant.gt_types)):
            #for gt in variant.gt_types:
                if(len(self.__my_target_indvs) > 0):
                    if(self.__my_target_indvs[i] == 1):
                        fields_list.append(self.__determine_gt_code(clas,variant.gt_types[i]))
                else:
                    fields_list.append(self.__determine_gt_code(clas, variant.gt_types[i]))
            #print "\n"
            line = " ".join(fields_list)

            #print "at the end of locus genotype file data line maker"
            #print line
            return line

    def __make_locus_name(self, prefix, counter,variant):
        new_name = "_".join((prefix,str(counter)))

        if(self.__LG_group_file != None):
            #do something to read data from LG_groups file and append it to the SNP name
            new_name = new_name + "_LG_" + self.__get_LG_group(variant)

        if(self.__subset_file != None):
            #do something to append frequency among the target individuals to the SNP name
            new_name = new_name + "_" + str(self.__get_subset_MAF(variant))

        return new_name

    def set_subset_file(self,filename):
        self.__subset_file = filename

    def set_LG_group_file(self, file):
        self.__LG_group_file = file

    def __get_LG_group(self,variant):

        LG_string = "no_LG"

        if(variant.CHROM.startswith("library") == True):
            return "R_gen"


        for line in open(self.__LG_group_file,'r'):

            line = line.strip()
            if(len(line) > 0 and line.startswith("#") == False
                    and line.startswith("MergeMark") == False):

                line_parts = line.strip().split(",")
                LG = line_parts[2]
                mapped_scaffold = line_parts[13]

                if(mapped_scaffold != "NA" and
                        variant.CHROM == mapped_scaffold):
                    LG_string = str(LG)


        return str(LG_string)

    def __get_subset_MAF(self,variant):

        gt_count = 0
        denominator = 0

        temporary_target_indvs = self.setup_group_member_list(self.__subset_file)

        for i in range(len(variant.gt_types)):
            if(temporary_target_indvs[i] == 1):
                if(variant.gt_types[i] == 3):
                    gt_count = gt_count + 2
                    denominator = denominator + 2
                elif(variant.gt_types[i] == 1):
                    gt_count = gt_count + 1
                    denominator = denominator + 2
                elif(variant.gt_types[i] == 0):
                    denominator = denominator + 2
        if(denominator > 0):
            #print "variant:\t" + variant.CHROM + "\t" + str(variant.start) + "\t" + str(variant.end)
            #print "frequency is:\t" + str(float(gt_count)/float(denominator))
            return round((float(gt_count) / float(denominator)),2)
        else:
            return 0

    def __determine_gt_code(self,clas_tuple,cyvcf2_gt_code):

        field_gt_code = ""
        if(len(clas_tuple) == 2 and self.__my_poptype == "HAP"):
            if(cyvcf2_gt_code == 0):
                field_gt_code = "b"
            elif(cyvcf2_gt_code == 2):
                field_gt_code = "-"
            elif(cyvcf2_gt_code == 3):
                field_gt_code = "a"
        elif(self.__my_poptype == "CP"):
            if (cyvcf2_gt_code == 0):
                field_gt_code = "kk"
            elif (cyvcf2_gt_code == 1):
                field_gt_code = "hk"
            elif (cyvcf2_gt_code == 2):
                field_gt_code = "--"
            elif (cyvcf2_gt_code == 3):
                field_gt_code = "hh"
            #if(len(clas_tuple) == 2):
            #    if(cyvcf2_gt_code == 0):
            #        field_gt_code = "ll"
            #    elif (cyvcf2_gt_code == 1):
            #        field_gt_code = "lm"
            #    elif(cyvcf2_gt_code == 2):
            #        field_gt_code = "--"
            #elif(len(clas_tuple) == 3):
            #    if(cyvcf2_gt_code == 0):
            #        field_gt_code = "kk"
            #    elif(cyvcf2_gt_code == 1):
            #        field_gt_code = "hk"
            #    elif(cyvcf2_gt_code == 2):
            #        field_gt_code = "--"
            #    elif(cyvcf2_gt_code == 3):
            #        field_gt_code = "hh"
        #print "cyvcf2_gt_code is:\t" + str(cyvcf2_gt_code)
        #print clas_tuple
        #print "field gt code is:\t" + str(field_gt_code)
        return field_gt_code

    def __get_hom_ref_code(self,clas_tuple):
            return clas_tuple[0]

    def __get_het_code(self,clas_tuple):
        return clas_tuple[1]

    def __get_hom_alt_code(self,clas_tuple):
        if(len(clas_tuple) == 4):
            return clas_tuple[3]
        else:
            return clas_tuple[1]

    def __clas_to_string(self,clas):

        clas_string = ",".join(clas)
        clas_string = "(" + clas_string + ")"
        return clas_string

    def __determine_segregation(self,gt_type_list):

        segregation=""
        clas=""
        phase=""

        het_count = 0
        hom_ref_count = 0
        hom_alt_count = 0
        missing_count = 0

        for i in range(len(gt_type_list)):
            if(len(self.__my_target_indvs) > 0):
                if(self.__my_target_indvs[i] ==1):
                    gt = gt_type_list[i]
                    if(gt == 0):
                        hom_ref_count = hom_ref_count + 1
                    elif(gt == 1):
                        het_count = het_count + 1
                    elif(gt == 2):
                        missing_count = missing_count + 1
                    elif(gt == 3):
                        hom_alt_count = hom_alt_count + 1
            else:
                gt = gt_type_list[i]
                if(gt == 0):
                    hom_ref_count = hom_ref_count + 1
                elif(gt == 1):
                    het_count = het_count + 1
                elif(gt == 2):
                    missing_count = missing_count + 1
                elif(gt == 3):
                    hom_alt_count = hom_alt_count + 1

        #print "het_count\thom_ref_count\thom_alt_count\tmissing_count"
        #print "\t".join((str(het_count),str(hom_ref_count),str(hom_alt_count),str(missing_count)))

        if(self.__my_poptype == "CP"):
            #if(hom_alt_count > 0 and het_count > 0 and hom_ref_count > 0):
            #    segregation = "<hkxhk>"
            #    clas = ["hh","hk","kk"]
            #    phase = "{11}"
            #elif(hom_alt_count > 0 or hom_ref_count > 0):
            #    segregation = "<lmxll>"
            #    clas = ["ll","lm"]
            #    phase = "{1-}"
            segregation = "<hkxhk>"
            clas = ["hh","hk","kk"]
            phase = "{11}"
            return segregation, clas, phase
        elif(self.__my_poptype == "HAP"):
            segregation = None
            clas = ["a","b"]
            phase = "{1}"
            return segregation, clas, phase

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",type=str,help="vcf file to prepare for joinmap")
    parser.add_argument("--name",type=str,help="name for the population. canno contain spaces")
    parser.add_argument("--popt",type=str,help="code for the population type. See Table 1 of JoinMap 4 maual\n/"
                                               "one of the following values: BC1, F2, RIx, DH, DH1, DH2, HAP\n/"
                        "HAP1, CP, BCpxFy, IMxFy. Case sensitive")
    parser.add_argument("--sample_list",type=str,help="name of a file with a list of samples to include")
    parser.add_argument("--interval_file",type=str,help="optional bed file that defines which intervals to use for the output file")
    parser.add_argument("--LG_map_file",type=str,help="optional file that contains a map to get linkage groups")
    parser.add_argument("--subset_samples",type=str,help="optional file for calculating MAFs ")
    parser.add_argument("--phenotype_file",type=str,help="name of a file with list of samples and phenotypes (binary, 1/0)")
    #add name for phenotype as a new argument
    #add param for methods to choose one SNP per scaffold
    args = parser.parse_args()
    print args
    main(args)

