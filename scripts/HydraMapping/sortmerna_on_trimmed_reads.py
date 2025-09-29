import subprocess
import os

'''sortmerna_process

    Function for executing sortmerna.

'''
def sortmerna_process(input_directory,sortme_output_dir, sample_name):
    # check if input and output paths exist
    if os.path.isdir(input_directory) == False or os.path.isdir(sortme_output_dir) == False:
        print("[-] ERROR either the specified input directory or the output directory (or both) do not exist!")
        return 1

    print("[*] Input and output directories exist")
    print("[*] Setting up logfile ...")
    logfile_path = sortme_output_dir + sample_name + "_sortmerna.log"

    with open(logfile_path,'w') as logfile:
        # check if rRNA databases exist
        logfile.write("INFO:checking if reference rRNA databases exist ...\n")
        reference_db_18s = "/gpfs/project/lubec100/CuTrOmics/rRNADBs/SILVA/silva-euk-18s-id95.fasta"
        reference_db_28s = "/gpfs/project/lubec100/CuTrOmics/rRNADBs/SILVA/silva-euk-28s-id98.fasta"
        if os.path.isfile(reference_db_18s) == False or os.path.isfile(reference_db_28s) == False:
            logfile.write("ERROR:reference database does not exist!\n")
            print("[-] ERROR either the specified input directory or the output directory (or both) do not exist!")
            return 1
        else:
            print("[+] REFERENCE databases exists ...")
            logfile.write("INFO:databses exist ...\n")

        logfile.write("INFO:extracting reads in {} directory ...\n".format(input_directory))

        # reading fastq files
        files = []
        for file in os.listdir(input_directory):
            if file.endswith(".fastq.gz") == True:
                files.append(file)

        # extracting sample names from files
        # creating dictionary with fw reads and rev reads as values for each sample (key)
        logfile.write("INFO:starting to extract samples with fw and rev reads ...\n")
        ordered_samples = {}
        samples = []

        for file in files:
            # the sample names reside in the file names
            if "fw_paired_out_" in str(file) or "r_paired_out_" in str(file):
                sample = file.split("_R")[0]

                logfile.write("\tINFO:found sample {}\n".format(sample))
                if "R1" in file:
                    sample = sample.split("fw_paired_out_")[1]
                elif "R2" in file:
                    sample = sample.split("r_paired_out_")[1]

                remainder = file.split("_R")[1]

                if sample not in samples:
                    samples.append(sample)
                    if remainder.startswith("1"):
                        reverse_remainder = "2" + remainder[1:]
                        ordered_samples[sample] = ["/fw_paired_out_"+sample+"_R"+remainder, "/r_paired_out_"+sample+"_R"+reverse_remainder]
                    elif remainder.startswith("2"):
                        forward_remainder = "1" + remainder[1:]
                        ordered_samples[sample] = ["/fw_paired_out_"+sample+"_R"+forward_remainder, "/r_paired_out_"+sample+"_R"+remainder]

        print("[*] length samples: {} --> total PE files to process {}".format(len(samples),len(files)))
        logfile.write("INFO:working with {} samples and {} PE files in total\n".format(len(samples),len(files)))

        #looping over samples to execute sortmerna
        for sample in samples:
            fw_file = input_directory + ordered_samples[sample][0]
            r_file = input_directory + ordered_samples[sample][1]

            if os.path.isfile(fw_file) == False or os.path.isfile(r_file) == False:
                print("[-] ERROR forward or reverse file not found for sample: {}".format(sample))
                logfile.write("ERROR:forward or reverse file not found for sample: {}\n".format(sample))
                return 1

            sample_output_rRNA = sortmerna_output_dir + sample + "_rRNA"
            sample_output_nonrRNA = sortme_output_dir + sample + "_nonrRNA"

            if os.path.isfile(sample_output_nonrRNA + "_fwd.fq.gz") == True and os.path.isfile(sample_output_nonrRNA + "_rev.fq.gz") == True:
                print("[*] Files for sample: {} already exist, skipping procedure ...".format(sample))
                logfile.write("INFO:Files for sample {} already exist, skipping ...\n".format(sample))
            else:
                # ensure to specify a working directory for every sample, sortmerna only works if the indexing directory kvdb is empty or does not exist!
                sortmerna_working_directory = sortme_output_dir + sample

                # example invocation
                '''
                sortmerna \
                  --ref ~/sortmerna_db/SSU_EUK \
                  --ref ~/sortmerna_db/LSU_EUK \
                  --reads sample_R1.fastq \
                  --reads sample_R2.fastq \
                  --aligned ~/sortmerna_results/sample_rRNA \
                  --other ~/sortmerna_results/sample_nonrRNA \
                  --fastx \
                  --paired_out \
                  --log ~/sortmerna_results/sortmerna.log
                '''

                number_threads = 8
                cmd = "sortmerna --ref {} --ref {} --reads {} --reads {} --aligned {} --other {} --fastx --threads {} --workdir {} --out2".format(reference_db_18s, reference_db_28s, fw_file, r_file,sample_output_rRNA, sample_output_nonrRNA,number_threads, sortmerna_working_directory)
                logfile.write("INFO:executing RSEM with following command: {}\n".format(cmd))
                proc = subprocess.Popen(cmd, shell=True)
                returncode = proc.wait()  # python 2 restriction, maybe better use run instead of Popen
                if(returncode != 0):
                    print("[-] ERROR during Popen cmd for sortmerna")
                    logfile.write("ERROR:Popen cmd for sortmerna error - returncode for Popen != 0 ...\n")
                    return 1
                else:
                    logfile.write("INFO:DONE sortmerna procedure for sample {}\n".format(sample))
                    print("[+] DONE sortmerna procedure for sample {}".format(sample))


                print("[*] INFO: removing temporary output ...\n")
                proc = subprocess.Popen("rm -r {}".format(sortmerna_working_directory), shell=True)

                returncode = proc.wait()
                if(returncode != 0):
                    print("[-] ERROR during Popen cmd for sortmerna")
                    logfile.write("ERROR:Popen cmd for rm error - returncode for Popen != 0 ...\n")
                    return 1
                else:
                    logfile.write("INFO:DONE removing directory for sample {}\n".format(sample))
                    print("[+] DONE removing directory for sample {}".format(sample))

        logfile.write("SUCESS\n")
    print("[+] DONE")
    return 0

sortmerna_output_dir = "/gpfs/project/lubec100/CuTrOmics/EcoKD1/stringent_trimmed_files/sortmerna_results/"
input_directory = "/gpfs/project/lubec100/CuTrOmics/EcoKD1/stringent_trimmed_files/"
sample_name = "EcoKD1"
sortmerna_process(input_directory,sortmerna_output_dir,sample_name)

sortmerna_output_dir = "/gpfs/project/lubec100/CuTrOmics/HydraAHL/stringent_trimmed_files/sortmerna_results/"
input_directory = "/gpfs/project/lubec100/CuTrOmics/HydraAHL/stringent_trimmed_files/"
sample_name = "HydraAHL"
sortmerna_process(input_directory,sortmerna_output_dir,sample_name)

sortmerna_output_dir = "/gpfs/project/lubec100/CuTrOmics/HydraRecolonization/trimmed_files/sortmerna_results/"
input_directory = "/gpfs/project/lubec100/CuTrOmics/HydraRecolonization/trimmed_files/"
sample_name = "HydraRecolonization"
sortmerna_process(input_directory,sortmerna_output_dir,sample_name)

sortmerna_output_dir = "/gpfs/project/lubec100/CuTrOmics/HydraTemperature/trimmed_files/sortmerna_results/"
input_directory = "/gpfs/project/lubec100/CuTrOmics/HydraTemperature/trimmed_files/"
sample_name = "HydraTemperature"
sortmerna_process(input_directory,sortmerna_output_dir,sample_name)
