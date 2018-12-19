import logging
import sys
import time
import requests
import os
import urllib.request 
import subprocess
import hashlib
import argparse

parser = argparse.ArgumentParser(description='Download a list of samples from ENA')
parser.add_argument('samplesheet', help='ENA sampelsheet containing accession ID, md5sums and download location.')
parser.add_argument('download_location', help='Location to store the downloaded fastq files at')
parser.add_argument('--aspera_binary', help='Location of the Aspera binary (default: use from PATH)', default='ascp')
parser.add_argument('--aspera_openssh', help='Location of the Aspera openssh (default: %(default)s)', default='~/.aspera/connect/etc/asperaweb_id_dsa.openssh')
parser.add_argument('--download_speed', help='Aspera download speed in MB/s (default: %(default)sMB/s)', default=2000)
parser.add_argument('--sample', help='Single sample to download. Overwrites inclusion and exclusion list')
parser.add_argument('--inclusion_list_file', help='Newline separated file with list of samples to include')
parser.add_argument('--exclusion_list_file', help='Newline separated file with list of samples to exclude')
parser.add_argument('-i','--include_all_samples', action='store_true', help='Raise error if sample is not found in samplesheet. If not set, give warning')


args = parser.parse_args()

# Original version from https://github.com/npklein/publicRNAseq

format = '%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format=format, datefmt="%Y-%m-%d %H:%M:%S")


class Download_ENA_samples:
    def __init__(self, samplesheet, download_location, aspera_binary='ascp', 
                 aspera_openssh='~/.aspera/connect/etc/asperaweb_id_dsa.openssh',
                 download_speed=2000,inclusion_list = [], exclusion_list = [],
                 include_all_samples = False):
        '''Initiate download_ENA_samples class
        
        samplesheet(str)    Samplesheet downloaded from http://www.ebi.ac.uk/ena/data/warehouse/search
                            Reports tab, with all columns selected.
        download_location(str):   Location to store the downloaded fastq files at
        aspera_binary(str)  Location of the Aspera binary (default: use from PATH) 
        aspera_openssh(str)  Location of the Aspera openssh (default: ~/.aspera/connect/etc/asperaweb_id_dsa.openssh) 
        inclusion_list(list)    Samples to download (def: [] -> no samples get excluded)
        exclusion_list(list)    Samples to exclude from download (def: None -> all samples get included)
        include_all_samples(Boolean) Check that all samples are in the samplesheet and downloaded, or throw error (default: False)
        '''
        # os.path.expanduser changes ~ into homedir
        self.aspera_binary = os.path.expanduser(aspera_binary)
        self.aspera_openssh = os.path.expanduser(aspera_openssh)
        if not os.path.isdir(download_location):
            logging.error('Download destination directory '+download_location+' does not exist')
            raise RuntimeError('Download destination directory '+download_location+' does not exist')
        self.download_location = download_location
        self.samplesheet = samplesheet
        self.inclusion_list = inclusion_list
        self.exclusion_list = exclusion_list
        self.aspera_download_speed = download_speed
        self.previous_percent = -1
        self.include_all_samples = include_all_samples        

    def __check_if_aspera_exists(self, aspera_binary):
        '''Check if aspera location given in aspera_binary exists or is in PATH
        
           aspera_binary(str):   Location of the Aspera binary
           aspera_openssh(str)  Location of the Aspera openssh
        '''
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
        fpath, fname = os.path.split(aspera_binary)
        if fpath:
            if is_exe(aspera_binary):
                return aspera_binary
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, aspera_binary)
                if is_exe(exe_file):
                    return exe_file
        logging.error(aspera_binary+' does not exist and is not found in PATH. Set correct location for binary and openssh key '+
                      'when initiating Download_ENA_samples class or with set function')        
        raise RuntimeError('Aspera binary not found')
    
    def set_aspera_binary(self, aspera_binary):
        '''Set path to aspera binary'''
        self.aspera_binary = aspera_binary
    
    def set_aspera_openssh(self, aspera_openssh):
        '''Set path to aspera openssh'''
        self.aspera_openssh = aspera_openssh   
        
    def set_inclusion_list(self, inclusion_list):
        '''Samples to include for download'''
        self.inclusion_list = inclusion_list
    
    def set_exclusion_list(self, exclusion_list):
        '''Samples to exclude for download'''
        self.exclusion_list = exclusion_list
    
    def set_aspera_speed(self, speed):
        '''Set aspera download speed
        
           speed(int)   Download speed for Aspera
        '''
        self.aspera_download_speed = speed
    
    def __reporthook(self,blocknum, blocksize, totalsize):
        ''' Report hook for urlretrieve, will print progress of download
        Copied from: http://stackoverflow.com/a/13895723/651779
        '''
        readsofar = blocknum * blocksize
        if totalsize > 0:
            percent = readsofar * 1e2 / totalsize
            s = "\r%5.1f%% %*d / %d" % (
                percent, len(str(totalsize)), readsofar, totalsize)
            if (round(percent) == 0 or round(percent) % 5 == 0) and self.previous_percent != round(percent):
                logging.info(s)
            if readsofar >= totalsize:
                if round(percent) % 5 == 0 and self.previous_percent != round(percent): 
                    logging.info("\n")
            if round(percent) % 5 == 0 and self.previous_percent != round(percent): 
                self.previous_percent = round(percent)
        else: # total size is unknown
            logging.info("read %d\n" % (readsofar,))
      
    def __download_sample_with_aspera(self, fastq_aspera_link, download_location):
        '''Download fastq file using aspera
        
           fastq_aspera_link(str)   Aspera link to download
           download_location(str)   Location to save downloaded file at
        '''
        # -T: disable encryption for maximum throughput
        # -l: max_rate (in this case in MB)
        # -i: private key file
        # End with download location
        command = [self.aspera_binary,'-QT','-P','33001', '-l'+str(self.aspera_download_speed)+'m', '-i', self.aspera_openssh,fastq_aspera_link,download_location]
        logging.info('Downloading with the terminal command:\n'+' '.join(command))
        subprocess.call(command)
        
    def __check_md5(self, fastq_file, md5, throw_error=True):
        '''Check if the md5sum of downloaded fastq file is correct
    
           fastq_file(str):   Location of fastq file to check md5 of
           md5(str):   md5sum to check if it is same as md5sum of fastq_file
           throw_error(Boolean): Throw error if fastq file does not exist or of md5 don't match
        '''
        if not os.path.exists(fastq_file):
            logging.info(fastq_file+' does not exist')
            if throw_error:
                raise RuntimeError('Fastq file '+fastq_file+' does not exist')
            return False
        logging.info('Calculating md5 of '+fastq_file)
        with open(fastq_file,'rb') as file_to_check:
            data = file_to_check.read()    
            md5_returned = hashlib.md5(data).hexdigest()
        logging.info('Calculated md5 = '+md5_returned+', comparing to '+md5)
        if md5_returned == md5:
            logging.info('md5 sums are the same')
        else:
            logging.info('fastq file md5sums are not the same. Should be '+md5+' but found '+md5_returned)
            if throw_error:
                raise RuntimeError('fastq file md5sums are not the same. Should be '+md5+' but found '+md5_returned)
        return md5_returned == md5
    
    def __report_number_of_fastq_files(self, run_accession, number_of_fastq_download_links):
        '''Simple print function that is used for aspera and ftp downloads'''
        if number_of_fastq_download_links == 1:
            logging.info(run_accession+' seems to be single end, downloading 1 fastq file')
        elif number_of_fastq_download_links == 2:
            logging.info(run_accession+' seems to be paired end, downloading 2 fastq file')
        elif number_of_fastq_download_links == 3:
            logging.warn(run_accession+' seems to be paired end, but has 3 download 2 files. Downloading all 3')
        else:
            logging.error('Expected 1-3 fastq files, got'+str(number_of_fastq_download_links)+'. Check your ena samplesheet')
            logging.error('Unexpected number of fastq files. Expected 1-3 fastq files, got'+str(number_of_fastq_download_links))    
        
    def __get_all_indices(self, list_to_index):
        '''Get the indexes of all the items in a list and put them in a dict with key: element, value: index
        
           list_to_index(list)    List to get index from all elements from
        '''
        list_indexes = {}
        i = 0
        for element in list_to_index:
            list_indexes[element] = i
            i += 1
        return list_indexes

    def start(self, download_protocol='aspera'):
        '''Download the samples using either aspera or ftp
        
           download_protocol(str):   Download protocol to use (def: aspera). Can only be aspera or ftp
        '''
        if download_protocol == 'aspera':
            self.__check_if_aspera_exists(self.aspera_binary)
            if not os.path.isfile(self.aspera_openssh):
                logging.error('openssh key for Aspera not found at '+self.aspera_openssh)
                raise RuntimeError('openssh key for Aspera not found at '+self.aspera_openssh)
            logging.info('Found aspera binary at '+self.aspera_binary)
        elif download_protocol != 'ftp':
            logging.error('download_protocol variable given to download_samples was '+download_protocol+', not aspera or ftp')
            raise RuntimeError('download protocol can only be aspera or ftp')
        included_samples = []
        print('Downloading samples to '+self.download_location)
        with open(self.samplesheet,'r', encoding='utf-8') as samplesheet_handle:
            samplesheet_header = samplesheet_handle.readline().split('\t')
            header_index = self.__get_all_indices(samplesheet_header)
            for line in samplesheet_handle:
                line = line.strip().split('\t')
                run_accession = line[header_index['run_accession']]
                # exclude list overrides include list
                # first check if self.inclusion_list is not empty, as include list should only be used when at least one sample is given
                if self.inclusion_list and run_accession in self.inclusion_list:
                    included_samples.append(run_accession)
                    if run_accession not in self.exclusion_list:
                        logging.info('Found '+run_accession+' to process')
                        if download_protocol == 'ftp':
                            raise RuntimeError('ftp not implemented yet')
                            fastq_ftp_links = line[header_index['fastq_ftp']].rstrip(';').split(';')
                            self.__report_number_of_fastq_files(run_accession, len(fastq_ftp_links))
                            # paired end data will have multiple download links
                            for fastq_ftp_link in fastq_ftp_links:
                                download_file_location = self.download_location+'/'+fastq_ftp_link.split('/')[-1]
                                x = 0
                                while not self.__check_md5(download_file_location, line[header_index['fastq_md5']], self.include_all_samples):
                                    fastq_ftp_link = 'ftp://'+fastq_ftp_link
                                    logging.info('Downloading '+fastq_ftp_link+' using ftp to '+download_file_location+'...')
                                    urllib.request.urlretrieve(fastq_ftp_link, download_file_location, self.__reporthook)
                                    x += 1
                                    if x == 10:
                                        self.logging.warning('Tried 10 times to download '+run_accession+' but md5sum never correct. Skipping')
                                    break
                        elif download_protocol == 'aspera':
                            fastq_aspera_links = line[header_index['fastq_aspera']].rstrip(';').split(';')
                            logging.info('Download from:'+'\n\t'.join(fastq_aspera_links))

                            self.__report_number_of_fastq_files(run_accession,len(fastq_aspera_links))
                            for index, fastq_aspera_link in enumerate(fastq_aspera_links):
                                download_file_location = self.download_location+'/'+fastq_aspera_link.split('/')[-1]
                                x = 0
#                                while not self.__check_md5(download_file_location, line[header_index['fastq_md5']], self.include_all_samples):
                                # Sometimes the download fails, so check if the fastq file is there. But remove it first in case this job has already run
                                if os.path.exists(download_file_location):
                                    os.remove(download_file_location)
                                while True:
                                    fastq_aspera_link = 'era-fasp@'+fastq_aspera_link
                                    logging.info('Downloading '+fastq_aspera_link+' to '+download_file_location+' using aspera...')
                                    self.__download_sample_with_aspera(fastq_aspera_link, download_file_location)
                                    if os.path.exists(download_file_location):
                                        break
                                    x += 1
                                    if x == 10:
                                        self.logging.warning('Tried 10 times to download '+run_accession+' but never connected.')
                                        raise RuntimeError('Tried 10 times to download '+run_accession+' but never connected.')
                                    sleep_time = 600
                                    print('Download failed, sleeping '+str(sleep_time)+' seconds then try again')
                                    time.sleep(sleep_time)
                                    
                                self.__check_md5(download_file_location, line[header_index['fastq_md5']].split(';')[index], self.include_all_samples)
                        else:
                            logging.error('Download protocol was not ftp or aspera')
                            raise RuntimeError('download protocol was not ftp or aspera')
        not_included_samples = [x for x in self.inclusion_list if x not in included_samples]
        if len(not_included_samples):
            if self.include_all_samples:
                raise RuntimeError('Not all samples from include list were present in the samplesheet '+self.samplesheet+'.\nMissing: '+'\t'.join(not_included_samples))
            else:
                logging.warn('Not all samples from include list were present in the samplesheet '+self.samplesheet+'.\nMissing: '+'\t'.join(not_included_samples))



if __name__ == "__main__":
    if args.sample:
        inclusion_list = [args.sample]
        exclusion_list = []
    else:
        if args.inclusion_list_file:
            with open(args.inclusion_list_file) as input_file:
                inclusion_list = input_file.read().split('\n')
        else:
            inclusion_list == []
        
        if args.exclusion_list_file:
            with open(args.exclusion_list_file) as input_file:
                exclusion_list = input_file.read().split('\n')
        else:
            exclusion_list == []
    download = Download_ENA_samples(args.samplesheet, args.download_location, args.aspera_binary, 
                                    args.aspera_openssh, args.download_speed, inclusion_list, exclusion_list,
                                    args.include_all_samples)

    download.start()
