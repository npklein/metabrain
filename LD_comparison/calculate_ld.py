import argparse
import requests

parser = argparse.ArgumentParser(description='Compare LD between metabrain and eqtlgen')
parser.add_argument('topSNPs', help='File with the top SNPs for metaBrain and eqtlGen')
parser.add_argument('ldlink_token', help='Token for API access ldlink (https://ldlink.nci.nih.gov/?tab=apiaccess)')

args = parser.parse_args()

with open(args.topSNPs) as input_file, open('metaBrain_eqtlGen_topSNP_LD.txt','w') as out:
    out.write("gene\ttopSNP_metaBrain\ttopSNP_eqtlGen\tD'\tR2\n")
    input_file.readline()
    x = 0
    for line in input_file:
        print(line)
        x += 1
        if x % 1000 == 0:
            print(x)
        line = line.strip()
        rs1 = line.split('\t')[1]
        rs2 = line.split('\t')[2]

        request_link = 'https://ldlink.nci.nih.gov/LDlinkRest/ldpair?var1='+rs1+'&var2='+rs2+'&pop=CEU%2BTSI%2BGBR%2BIBS&token='+args.ldlink_token
        ldlink_result = requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldpair?var1='+rs1+'&var2='+rs2+'&pop=CEU%2BTSI%2BGBR%2BIBS&token='+args.ldlink_token)
        try:
            d = ldlink_result.text.split("D': ")[1].split("\n")[0]
            r = ldlink_result.text.split("R2: ")[1].split("\n")[0]
        except:
            if 'error":' in ldlink_result.text:
                d = 'error: '+ldlink_result.text.split('error": ')[1].split('"')[1]
                r = 'error: '+ldlink_result.text.split('error": ')[1].split('"')[1]
            else:
                print(request_link)
                print(rs1, rs2)
                print(ldlink_result.text)

                raise

        out.write(line.split('\t')[0]+'\t'+rs1+'\t'+rs2+'\t'+d+'\t'+r+'\n')
