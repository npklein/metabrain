import argparse
import gzip
parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('eQTL_file',
                    help='File with eQTL data')
parser.add_argument('gwas_file',
                    help='File with gwas data')
parser.add_argument('gwas_annotation_file',
                    help='File with gwas annotation')
parser.add_argument('tested_snps',
                    help='File with tested SNPs')
parser.add_argument('outfile',
                    help='outfile name')
args = parser.parse_args()

def main():
    SNPs_tested = read_tested_SNPs()
    trait_info = read_gwas_trait_info()
    traits, SNPs_tested_in_gwas = parse_gwas(SNPs_tested)
    parse_eQTL_and_write_output(traits, SNPs_tested_in_gwas)

def read_tested_SNPs():
    print('start reading tested SNPs')
    with gzip.open(args.tested_snps,'rt') as input_file:
        return(set(input_file.read().split('\n')))

def read_gwas_trait_info():
    trait_info = {}
    print('start reading GWAS annotation file')
    with open(args.gwas_annotation_file) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            trait_info[line[0].lower()] = [line[1], line[2]]
    return(trait_info)

def parse_sample_size_from_text(sample_size_text):
    # search for all numbers in the SAMPLE_SIZE column and sum them (e.g. 1434 cases 431 controls -> 1865)
    # the replace of "cases," to "cases, " helps incase it is like "cases,1934", probably should need more generic function as same might happen with other words
    word_list = sample_size_text.replace('cases,', 'cases, ').replace('controls,', 'controls, ').replace(',','').split()
    sample_size = 0
    # add all the values if they are digit and if not the words around is from ... families or type ... diabetes
    # and fix some other problems
    contains_twin_pair = False
    for index, word in enumerate(word_list):
        if word.isdigit():
            if index > 0:
                if not word_list[index-1] == 'from' and not word_list[index-1].lower() == 'type':
                    sample_size += int(word)
            else:
                sample_size += int(word)
        if index < len(word_list)-1:
            if 'twin' in word and 'pair' in word_list[index+1]:
                contains_twin_pair = True
    if contains_twin_pair:
        sample_size *= 2
    return(sample_size)

def parse_confidence_interval(confidence_interval_text):
    # if value is Odds Ratio or Beta value by looking at characters after unit range of confidence interval,
    #e.g. [1.23-1.54] ml mediconfidence_intervalne -> Beta, [1.23-1.54] -> OR
    if '[' not in confidence_interval_text:
        confidence_interval = confidence_interval_text.split(']')[0]
    if '[' not in confidence_interval_text:
        
    elif ']' not in confidence_interval_text:
        confidence_interval = confidence_interval_text.split('[')[1]
    else:
        confidence_interval = confidence_interval_text.split('[')[1].split(']')[0]

    if confidence_interval == 'NR':
        return 'NoType', str(OR_or_BETA), 'NoZscore'

    if ']' in confidence_interval_text:
        if len(confidence_interval_text.split(']')[1].strip()) > 2:
            beta_OR_type = 'BETA'
        else:
            beta_OR_type = 'OR'
    elif ' ' in confidence_interval_text.split:
        if len(confidence_interval_text.split(' ')[1].strip()) > 2:
            beta_OR_type = 'BETA'
        else:
            beta_OR_type = 'OR'
    else:
        beta_OR_type = 'OR'

    # confidence interval is of pattern \[\d+-\d+\].*    extract lower and upper bound
    # can be negative in confidence interval, so can't just split
    if confidence_interval == 'NR':
        return 'NoType', str(OR_or_BETA), 'NoZscore'
    confidence_interval_sign = confidence_interval[0]
    confidence_interval = confidence_interval.lstrip('(').rstrip(')')
    # split only on first incase there is a negative [-0.3--0.1]
    # some have weird minus signs in them (like \xe2\x80\x90), check on them and use for split
    if '\xe2\x80\x93' in confidence_interval:
        confidence_interval = confidence_interval.lstrip('-').split('\xe2\x80\x93',1)
    elif '\xe2\x80\x90' in confidence_interval:
        confidence_interval = confidence_interval.lstrip('-').split('\xe2\x80\x90',1)
    else:
        confidence_interval = confidence_interval.lstrip('-').split('-',1)

    if len(confidence_interval) == 1:
        confidence_interval = confidence_interval[0].split(',')

    confidence_interval_lower = float(confidence_interval[0])
    # clearning up confidence interval
    confidence_interval_upper = float(confidence_interval[1].lstrip('.').replace(',','.'))
    if confidence_interval_sign == '-':
        confidence_interval_lower = -confidence_interval_lower

    return beta_OR_type, confidence_interval_lower, confidence_interval_upper

def convert_beta_or_OR_to_zscore(OR_or_BETA, confidence_interval_text, sample_size):
    if len(OR_or_BETA.strip()) == 0:
        return 'NoType', str(OR_or_BETA), 'NoZscore'
    try:
        OR_or_BETA = float(OR_or_BETA)
    except ValueError:
        print(OR_or_BETA)
        raise

    beta_OR_type, confidence_interval_lower, confidence_interval_upper = parse_confidence_interval(confidence_interval_text)

    # sometimes the `OR or BETA` value is same as lower or upper confidence interval, which calculates to infinite Z-score. If this happens, change to average of interval
    if OR_or_BETA == confidence_interval_lower or OR_or_BETA == confidence_interval_upper:
        OR_or_BETA = (confidence_interval_lower+confidence_interval_upper)/2

    #calculate Z-score
    #taken from http://stats.stackexchange.com/a/69403 and https://www.biostars.org/p/140292/
    # the difference between BETA and confidence_interval_LOWER is needed for calculating the standard error (se)
    beta_conf_diff = OR_or_BETA - confidence_interval_lower
    # becuase it can happen dat the standard error is 0 which makes it impossible to calc. zscore
    # add a little bit to standard error to make sure that does not happen
    standard_error = (beta_conf_diff/1.96)+0.0000000001
    zscore = OR_or_BETA / standard_error
    zscore_weighted = zscore/float(sample_size)
    return(str(beta_OR_type), str(OR_or_BETA), str(zscore))

def parse_gwas(SNPs_tested):
    traits = {}
    SNPs_in_gwas = set([])
    print('start reading GWAS file')
    snp_gwas_info = {}
    with open(args.gwas_file) as input_file:
        header = input_file.readline().strip().split('\t')
        snps_index = header.index('SNPS')
        pval_index = header.index('P-VALUE')
        OR_or_beta_index = header.index('OR or BETA')
        risk_allele_index = header.index('STRONGEST SNP-RISK ALLELE')
        print(header)
        conf_interval_index = header.index('95% CI (TEXT)')
        sample_size_index = header.index('INITIAL SAMPLE SIZE')
        for line in input_file:
            line = line.strip().split('\t')
            trait = line[7]
            pval = line[pval_index]
            OR_or_beta = line[OR_or_beta_index]
            risk_allele = line[risk_allele_index].split('; ')
            conf_interval = line[conf_interval_index]
            sample_size = parse_sample_size_from_text(line[sample_size_index])
            SNPs = line[snps_index].split(';')
            try:
                beta_OR_type, OR_or_BETA, zscore = convert_beta_or_OR_to_zscore(OR_or_beta, conf_interval, sample_size)
            except:
                print(line)
                raise
            if not zscore:
                zscore = 'NoZscore'
            for index, snp in enumerate(SNPs):
                if 'rs' not in snp:
                    continue
                if snp not in traits:
                    traits[snp] = []
                traits[snp].append([trait, pval, OR_or_beta, risk_allele[index], conf_interval, beta_OR_type,  zscore])
                SNPs_in_gwas.add(snp)

    SNPs_tested_in_gwas = SNPs_tested & SNPs_in_gwas
    return(traits, SNPs_tested_in_gwas)


def parse_eQTL_and_write_output(traits, SNPs_tested_in_gwas):
    eqtl_SNPs = {}
    snp_seen = set([])
    highest_pvalue = -1
    print('start reading eQTL file '+args.eQTL_file)
    print('start writing results to '+args.outfile)
    with gzip.open(args.eQTL_file,'rt') as input_file, open(args.outfile,'wt') as out:
        out.write('SNP\tlowest_pval\ttrait\tparent_term\tmodified_parent_term\tgwas_pval\tgwas_OR_or_beta\tgwas_risk_allele\teqtl_allele_assessed\teqtl_zscore\tgwas_conf_interval\tgwas_OR_or_beta_type\tgwas_zScore\n')
        header = input_file.readline().strip('\n')
        split_head = header.split('\t')
        pval_index = split_head.index('PValue')
        snp_index = split_head.index('SNPName')
        allele_assessed_index = split_head.index('AlleleAssessed')
        zscore_index = split_head.index('OverallZScore')
        x = 0
        y = 0
        for line in input_file:
            x += 1
            if x == 100000:
                print(x,'lines read,',y,'lines written')
            line = line.strip('\n')
            split_line = line.split('\t')
            snp = split_line[snp_index]
            pval = split_line[pval_index]
            allele_assessed = split_line[allele_assessed_index]
            zscore = split_line[allele_assessed_index]
            highest_pvalue = pval
            if 'rs' in snp and ':' in snp:
                snp = snp.split(':')[2]
            if snp in snp_seen: # or snp not in SNPs_tested_in_gwas:
                continue
            snp_seen.add(snp)
            if snp not in traits:
                out.write(snp+'\t'+pval+'\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\t'+allele_assessed+'\t'+zscore+'NotInGWAS\tNotInGWAS\tNotInGWAS\n')
                continue
            traits_of_snp = traits[snp]
            for trait_gwas_info in traits_of_snp:
                trait = trait_gwas_info[0].lower()
                try:
                    parent_term = trait_info[trait][0]
                    modified_parent_term = trait_info[trait][1]
                except KeyError:
                    print(trait,'not found')
                    parent_term = 'UNKNOWN'
                    modified_parent_term = 'UNKNOWN'
                out.write(snp+'\t'+pval+'\t'+trait_gwas_info[0]+'\t'+parent_term+'\t'+modified_parent_term+'\t'+trait_gwas_info[1]+'\t')
                out.write(trait_gwas_info[2]+'\t'+trait_gwas_info[3]+'\t'+allele_assessed+'\t'+zscore+'\t'+trait_gwas_info[4]+'\t')
                out.write(traits_gwas_info[5]+'\t'+traits_gwas_info[6]+'\n')
                y += 1
        print('following list of SNPs not in eqtl file, giving the highest p-value ('+str(highest_pvalue)+')')
        for snp in SNPs_tested_in_gwas - snp_seen:
            if snp not in traits:
                out.write(snp+'\t'+pval+'\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tUNKNOWN\tUNKNOWN\tNotInGWAS\tNotInGWAS\tNotInGWAS\n')
                continue
            traits_of_snp = traits[snp]
            for trait_gwas_info in traits_of_snp:
                trait = trait_gwas_info[0].lower()
                try:
                    parent_term = trait_info[trait][0]
                    modified_parent_term = trait_info[trait][1]
                except KeyError:
                    print(trait,'not found')
                    parent_term = 'UNKNOWN'
                    modified_parent_term = 'UNKNOWN'
                out.write(snp+'\t'+pval+'\t'+trait_gwas_info[0]+'\t'+parent_term+'\t'+modified_parent_term+'\t'+trait_gwas_info[1])
                out.write('\t'+trait_gwas_info[2]+'\t'+trait_gwas_info[3]+'\tUNKNOWN\tUNKNOWN\t'+trait_gwas_info[4])
                out.write('\t'+trait_gwas_info[5]+'\t'+traits_gwas_info[6]+'\n')
        
    print('written to '+args.outfile)
    print('Done!')

if __name__ == "__main__":
    main()


