import argparse
import re

def main(args):
	header,table = parse_table(args.table)
	table = append_vep_to_table(table,args.vep)
	table = append_gff_alias_to_table(table, args.gff3)

	print_table_with_header(header,table)

def print_table_with_header(header,table):

	print header + "\timpact\ttranscriptID"
	for gene in table.keys():
		line = str(table[gene]["table_line"]) + "\t" + str(table[gene]["vep_impact"]) + "\t" + str(table[gene]["gene_ID"]) + "\t" + str(table[gene]["gff_alias"])
		print line


def append_gff_alias_to_table(table, gff3_file):
	gff3_file = open(gff3_file,'r')

	for line in gff3_file:
		m = re.search("gene",line)
		if(m):
			parts = line.strip().split("=")
			alias = parts[len(parts) - 1]
			ID = line.strip().split("\t")[8].split("=")[1].split(";")[0]

			for key in table.keys():
				if table[key]["gene_ID"] == ID:
							table[key]["gff_alias"] = alias
	return table

def parse_table(table_file):
	first_line = ""
	table = open(table_file,'r')

	SNP_table = {}
	for line in table:
		if(len(first_line) == 0):
			first_line = line.strip()
		else:
			line_parts = line.strip().split()
			table_map = {"table_line":0,"vep_impact":0,"gff_alias":0,"gene_ID":0}
			table_map["table_line"] = line.strip()
			SNP_table.setdefault(line_parts[1],table_map)
	table.close()
	return first_line,SNP_table

def append_vep_to_table(table,vep_file):

	vep = open(vep_file,'r')
	for line in vep:
		line_parts = line.strip().split()
		if(line_parts[1] in table):
			table[line_parts[1]]["vep_impact"] = line_parts[6].strip()
			table[line_parts[1]]["gene_ID"] = line_parts[3]
			#print table[line_parts[1]] + "\t" + line_parts[6].strip() + "\t" + line_parts[4].strip()
	vep.close()
	return table

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--table",type=str,help="file with the association data")
	parser.add_argument("--vep",type=str,help="vep.txt file")
	parser.add_argument("--gff3",type=str,help="gff3 file")
	args = parser.parse_args()
	main(args)
