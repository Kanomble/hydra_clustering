import subprocess
import os


'''collect_sortmerna_output

    Function to collect all output nonrRNA files of sortmerna for one particular dataset.
    The function outputs a dictionary with sample to fw,rev fastq file mapping.
    This mapping already contains the correct path, thus no adjustment is needed for the path.

'''
def collect_sortmerna_output(sortmerna_workdir):
    try:
        # reading fastq files
        files = []
        for file in os.listdir(sortmerna_workdir):
            if file.endswith(".fq.gz") == True:
                files.append(file)

        # extracting sample names from files
        sortmerna_sample_to_files = {}
        samples = []
        for file in files:
            # the sample names reside in the file names
            # by splitting the filenames we will get all
            # correct samples
            if "nonrRNA" in str(file):
                sample = file.split("_nonrRNA_")[0]
                if sample not in list(sortmerna_sample_to_files.keys()):
                    samples.append(sample)
                    fw_read = sortmerna_workdir+sample+"_nonrRNA_fwd.fq.gz"
                    rev_read = sortmerna_workdir+sample+"_nonrRNA_rev.fq.gz"
                    # check if files truly exist if not raise an error
                    if os.path.isfile(fw_read) == False or os.path.isfile(rev_read) == False:
                        raise Exception("Forward ({}) or reverse ({}) SortMeRNA output does not exist!".format(fw_read,rev_read))
                    else:
                        sortmerna_sample_to_files[sample] = [fw_read,rev_read]
        return sortmerna_sample_to_files
    except Exception as e:
        print("[-] ERROR during SortMeRNA sample to file mapping with exception: {}".format(e))
        raise Exception("[-] ERROR during SortMeRNA sample to file mapping with exception: {}".format(e))

'''repair_procedure

    SortMeRNA is introducing malformatted fastq files if one read aligns to rRNA DBs and the paired
    read does not align. Using the bbmap repair.sh program will fix this problem. Singleton output
    files can be processed later. The function takes a dictionary with a sample to files mapping as input.

'''
def repair_procedure(sortmerna_sample_to_files,repair_output_directory):
    try:
        # starting to fix files with bbmap/repair.sh
        if os.path.isdir(repair_output_directory) == False:
            print("[-] ERROR repair output directory does not exist")
            raise Exception("The repair output directory: {} does not exist!".format(repair_output_directory))

        repair_sample_to_files = {}
        for sample in sortmerna_sample_to_files.keys():
            print("[*] Repairing sample: {}".format(sample))
            fw_output = repair_output_directory + sample + "_repaired_fw.fq.gz"
            rev_output = repair_output_directory + sample + "_repaired_rev.fq.gz"
            singleton_output = repair_output_directory + sample + "_singletons.fastq"
            repair_sample_to_files[sample] = [fw_output,rev_output]
            if os.path.isfile(fw_output) == True and os.path.isfile(rev_output) == True:
                print("[*] Reformat output files already exist, skipping repair".format(sample))
            else:
                fw_input = sortmerna_sample_to_files[sample][0]
                rev_input = sortmerna_sample_to_files[sample][1]
                cmd = "~/bbmap/repair.sh in1={} in2={} out1={} out2={} outs={} -Xmx32g".format(fw_input,rev_input,fw_output,rev_output,singleton_output)
                print("[*] Starting reformatting for sample {} with cmd: {}".format(sample, cmd))
                proc = subprocess.Popen(cmd, shell=True)
                returncode = proc.wait()
                if(returncode != 0):
                    print("[-] ERROR during Popen cmd for REFORMAT")
                    raise Exception("ERROR during Popen cmd for repair.sh returncode != 0")
                else:
                    print("[*] DONE repairing sample: {}".format(sample))
        return repair_sample_to_files
    except Exception as e:
        print("[-] ERROR during repair procedure with exception: {}".format(e))
        raise Exception("[-] ERROR during repair procedure with exception: {}".format(e))

'''star_rsem_procedure

    This function performs the rsem-calculate-expression command. It takes an dictionary with
    a sample to file mapping, the files should be paired end files with the fw read as first
    argument and the rev read as second.

'''
def star_rsem_procedure(repair_sample_to_files,rsem_output_directory,star_genome):
    try:
        # looping over samples
        rsem_result_dict = {}
        for sample in repair_sample_to_files.keys():

            if os.path.isfile(rsem_output_directory + sample + '.genes.results') == False:
                print("[*] Output files do not exist, starting RSEM procedure for sample: {}".format(sample))
                fw_file = repair_sample_to_files[sample][0]
                r_file = repair_sample_to_files[sample][1]
                if os.path.isfile(fw_file) == False or os.path.isfile(r_file) == False:
                    print("[-] ERROR forward or reverse file not found for sample: {}".format(sample))
                    raise Exception("ERROR during star_rsem_procedure, fw_file {} or rev_file {} not exist!".format(fw_file, r_file))
                else:
                    # RSEM procedure
                    rsem_output = rsem_output_directory + sample
                    print("[*] Starting RSEM procedure for sample : {} ...".format(sample))
                    # rsem-calculate-expression -p 8 --paired-end --star --star-gzipped-read-file /gpfs/project/lubec100/CuTrOmics/EcoKD1/stringent_trimmed_files/deduplicated_reads/dedup_forward_read_I25740-L1_S30_L002.fastq.gz /gpfs/project/lubec100/CuTrOmics/EcoKD1/stringent_trimmed_files/deduplicated_reads/dedup_reverse_read_I25740-L1_S30_L002.fastq.gz /gpfs/project/lubec100/CuTrOmics/RSEMResults/index_files/HVAEP_reference
                    # /gpfs/project/lubec100/CuTrOmics/RSEMResults/expression_results
                    # constructing STAR cmd for Popen
                    number_threads = 8
                    cmd = "rsem-calculate-expression -p {} --paired-end --star --star-gzipped-read-file {} {} {} {}".format(number_threads, fw_file, r_file, star_genome,  rsem_output)
                    proc = subprocess.Popen(cmd, shell=True)
                    returncode = proc.wait()  # python 2 restriction, maybe better use run instead of Popen
                    if(returncode != 0):
                        print("[-] ERROR during Popen cmd for rsem-calculate-expression")
                        raise Exception("ERROR during Popen cmd for rsem-calculate-expression")
                    else:
                        rsem_result_dict[sample] = rsem_output
                        print("[+] DONE rsem-calculate-expression procedure of sample {}".format(sample))
            else:
                print("[+] FILE: {} already exist, skipping rsem-calculate-expression task!".format(rsem_output_directory + sample + '.genes.results'))
        return rsem_result_dict
    except Exception as e:
        print("[-] ERROR during rsem-calculate-expression taks with exception: {}".format(e))
        raise Exception("[-] ERROR during rsem-calculate-expression task with exception: {}".format(e))

'''main_repair_star_rsem_procedure

    Function for executing RSEM with the STAR aligner.

'''
def main_repair_star_rsem_procedure(sortmerna_input_directory,repair_output_directory,rsem_output_directory, sample_name):
    print("[*] Starting STAR alignment procedure with sample {}".format(sample_name))
    print("[*] Producing output files in directory : {}".format(rsem_output_directory))

    if os.path.isdir(sortmerna_input_directory) == False or os.path.isdir(rsem_output_directory) == False or os.path.isdir(repair_output_directory) == False:
        print("[-] ERROR either the specified input directory or the output directory (or both) do not exist!")
        raise Exception("ERROR one of the specified input directories does not exist!")

    print("[*] Input and output directories exist")
    print("[*] Setting up logfile ...")
    logfile_path = rsem_output_directory + sample_name + "_rsem_aligning.log"

    with open(logfile_path,'w') as logfile:
        try:
            # check if STAR genome index exist
            genome_name = "HVAEP_reference"
            star_genome_dir = "/gpfs/project/lubec100/CuTrOmics/RSEMResults/index_files/"
            if os.path.isdir(star_genome_dir) == False:
                print("[-] STAR genome directory : {} does not exist".format(star_genome_dir))
                raise Exception("ERROR STAR genome directory {} does not exist!".format(star_genome_dir))

            # collect sortmerna output for the repair program
            logfile.write("INFO:starting file collection\n")
            sortmerna_sample_to_files_dict = collect_sortmerna_output(sortmerna_input_directory)
            # perform repair and collect output files
            logfile.write("INFO:starting repair procedure with {} samples\n".format(len(sortmerna_sample_to_files_dict.keys())))
            repair_sample_to_files_dict = repair_procedure(sortmerna_sample_to_files_dict,repair_output_directory)
            logfile.write("INFO:starting rsem-calculate-expression procedure with {} samples\n".format(len(repair_sample_to_files_dict)))
            # perform star-rsem mapping procedure
            star_genome = star_genome_dir + genome_name
            rsem_output_dict = star_rsem_procedure(repair_sample_to_files_dict,rsem_output_directory,star_genome)
            logfile.write("INFO:DONE rsem-calculate-expression\n")
            return 0
        except Exception as e:
            logfile.write("ERROR:error during STAR-RSEM procedure with exception: {}".format(e))
            print("[-] ERROR during STAR-RSEM procedure with exception: {}".format(e))
            raise Exception("[-] ERROR during STAR-RSEM procedure with exception: {}".format(e))

# execution of functions for the different hydra bulk seq transcriptomics

rsem_output_dir = "/gpfs/project/lubec100/CuTrOmics/RSEMResultsAfterSortMeRNA/EcoKD1/"
repair_output_dir = "/gpfs/project/lubec100/CuTrOmics/EcoKD1/stringent_trimmed_files/sortmerna_results/cleaned_reads/"
sortmerna_input_dir = "/gpfs/project/lubec100/CuTrOmics/EcoKD1/stringent_trimmed_files/sortmerna_results/"
sample_name = "EcoKD1"
main_repair_star_rsem_procedure(sortmerna_input_dir,repair_output_dir,rsem_output_dir,sample_name)

rsem_output_dir = "/gpfs/project/lubec100/CuTrOmics/RSEMResultsAfterSortMeRNA/HydraAHL/"
repair_output_dir = "/gpfs/project/lubec100/CuTrOmics/HydraAHL/stringent_trimmed_files/sortmerna_results/cleaned_reads/"
sortmerna_input_dir = "/gpfs/project/lubec100/CuTrOmics/HydraAHL/stringent_trimmed_files/sortmerna_results/"
sample_name = "HydraAHL"
main_repair_star_rsem_procedure(sortmerna_input_dir,repair_output_dir,rsem_output_dir,sample_name)

rsem_output_dir = "/gpfs/project/lubec100/CuTrOmics/RSEMResultsAfterSortMeRNA/HydraRecolonization/"
repair_output_dir = "/gpfs/project/lubec100/CuTrOmics/HydraRecolonization/trimmed_files/sortmerna_results/cleaned_reads/"
sortmerna_input_dir = "/gpfs/project/lubec100/CuTrOmics/HydraRecolonization/trimmed_files/sortmerna_results/"
sample_name = "HydraRecolonization"
main_repair_star_rsem_procedure(sortmerna_input_dir,repair_output_dir,rsem_output_dir,sample_name)

rsem_output_dir = "/gpfs/project/lubec100/CuTrOmics/RSEMResultsAfterSortMeRNA/HydraTemperature/"
repair_output_dir = "/gpfs/project/lubec100/CuTrOmics/HydraTemperature/trimmed_files/sortmerna_results/cleaned_reads/"
sortmerna_input_dir = "/gpfs/project/lubec100/CuTrOmics/HydraTemperature/trimmed_files/sortmerna_results/"
sample_name = "HydraTemperature"
main_repair_star_rsem_procedure(sortmerna_input_dir,repair_output_dir,rsem_output_dir,sample_name)
