import argparse
from subprocess import run
from cyvcf2 import VCF, Writer
import os
import sys

def alt_present(pileup_output_path, alt):
    with open(pileup_output_path, 'r') as f:
        for line in f:
            var_piles = line.rstrip().split()[4]
            # If the alternate allele in present in even one of the reads
            if alt in var_piles:
                return True
            else:
                return False


def evaluate_dnv(dnv_info, ref_fa_path, mom_bam_path, dad_bam_path):
    dump = open('/dev/null', 'w')
    region = '{}:{}-{}'.format(dnv_info['chrom'], dnv_info['pos'], dnv_info['pos'])

    # samtools pileup at variant position in mother's BAM
    pileup_args = ["samtools", "mpileup", "-r", region, "-f", ref_fa_path, mom_bam_path]
    mom_pileup_ouput_path = os.path.join("mom.txt")
    with open(mom_pileup_ouput_path, 'w') as file:
        run(pileup_args, check=True, stdout=file,
            stderr=dump, universal_newlines=True)
        
    if alt_present(mom_pileup_ouput_path, dnv_info['alt'][0]):
        return "FP"
    
    # samtools pileup at variant position in father's BAM
    else:       
        pileup_args = ["samtools", "mpileup", "-r", region, "-f", ref_fa_path, dad_bam_path]
        dad_pileup_ouput_path = os.path.join("dad.txt")
        with open(dad_pileup_ouput_path, 'w') as file:
            run(pileup_args, check=True, stdout=file,
                stderr=sys.stderr, universal_newlines=True)
        if not alt_present(dad_pileup_ouput_path, dnv_info['alt'][0]):
            return "TP"
        else:
            return "FP"   


parser = argparse.ArgumentParser()
parser.add_argument('-mbam', help='BAM file (indexed)', required=True)
parser.add_argument('-dbam', help='BAM file (indexed)', required=True)
parser.add_argument('-r', help='reference FASTA file (indexed)', required=True)
parser.add_argument('-v', help='VCF file (can be bgzipped)',required=True)
parser.add_argument('-o', default='out.vcf', help='output (annotated) VCF (will be bgzipped if ending in .gz)')
# parser.add_argument('-d', help='output directory', default='temp_valsv')
args = parser.parse_args()

vcf = VCF(args.v)
vcf.add_info_to_header({'ID': 'flag', 
                        'Description': 'Flagging True-positive denovos identifed with in-silico denovo validation',
                        'Type': 'Float', 'Number': '1'})

vcf_o = Writer(args.o, vcf)

# Read VCF and evaluate each denovo SNP
for variant in vcf:
    dnvinfo = {}
    dnvinfo['chrom'] = variant.CHROM
    dnvinfo['pos'] = variant.POS
    dnvinfo['alt'] = variant.ALT

    flag = evaluate_dnv(dnvinfo, args.r, args.mbam, args.dbam)
    if flag == 'TP':
        variant.INFO['flag'] = 'TP'

    vcf_o.write_record(variant)


# os.path.join(output_dir, sv_info['svid'] + ".vcf")
bgzip_args = ["bgzip", "-c", vcf_o]
vcf_path_gz = vcf_o + '.gz'

vcf_o.close()
vcf.close()