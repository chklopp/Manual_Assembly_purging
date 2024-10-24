import yaml
import os
import sys
import argparse
import shutil
import subprocess
import numpy as np
from collections import defaultdict
import pickle


# Function to load YAML file into a Python dictionary
def load_yaml_config(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

# Function to check if a file exists on the file system
def is_file_accessible(file_path):
    return os.path.exists(file_path)

# Function to check the 'reads' and 'assemblies' files
def check_files(config):
    # Check the 'reads' file
    reads_file_path = config.get('reads')
    if reads_file_path:
        if is_file_accessible(reads_file_path):
            print(f"The reads file '{reads_file_path}' is accessible.")
        else:
            print(f"The reads file '{reads_file_path}' is not accessible or does not exist.", file=sys.stderr)
            sys.exit(1)  # Stop execution on failure

    # Check each file in the 'assemblies' list
    assemblies = config.get('assemblies', [])
    for assembly_file in assemblies:
        if is_file_accessible(assembly_file):
            print(f"The assembly file '{assembly_file}' is accessible.")
        else:
            print(f"The assembly file '{assembly_file}' is not accessible or does not exist.", file=sys.stderr)
            sys.exit(1)  # Stop execution on failure

# Function to parse a FASTA file and return a list of tuples (sequence_name, sequence_length)
def parse_fasta(assembly_file):
    sequences = []
    with open(assembly_file, 'r') as file:
        sequence_name = None
        sequence_length = 0
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                # If a sequence was being processed, store it before moving to the next one
                if sequence_name:
                    sequences.append((sequence_name, sequence_length))
                sequence_name = line[1:].split()[0]  # Take the first word as the sequence name
                sequence_length = 0  # Reset the sequence length for the new sequence
            else:
                # Add the length of the sequence line to the current sequence
                sequence_length += len(line)
        # Don't forget to add the last sequence after the loop
        if sequence_name:
            sequences.append((sequence_name, sequence_length))
    
    return sequences

# Function to parse a FASTA file like count file and return an array with sequence names and content
def parse_count(count_file):
    counts = []
    with open(count_file, 'r') as file:
        sequence_name = None
        sequence = "" 
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                # If a sequence was being processed, store it before moving to the next one
                if sequence_name:
                    counts.append((sequence_name, sequence.split(" ")))
                sequence_name = line[1:].split()[0]  # Take the first word as the sequence name
                sequence = ""  # Reset the sequence length for the new sequence
            else:
                # Add the length of the sequence line to the current sequence
                sequence += line
        # Don't forget to add the last sequence after the loop
        if sequence_name:
            counts.append((sequence_name, sequence.split(" ")))
    #print(counts)
    return counts

# Function to count the number of kmer of a certain class in the assembly file 
def count_per_slice(loc, sequences, config, slices):
    slice_size = int(config.get('slice_size'))  
    loc = loc.replace(".locate","")
    #loc = loc.rsplit('.', maxsplit=1)[0]
    
    # for each contig
    for (name, seq) in sequences:
        # create the empty numpy table which will recieve the count
        numpy_table = np.zeros((int(len(seq)/slice_size)+1, len(config.get('assemblies'))))

        # for each chunck and each kmer occurrences count de presence of this kmer count
        for i in range(0, int(len(seq)/slice_size)+1):
            chunk = seq[i*slice_size:i*slice_size+slice_size-1]
            for a in range(len(config.get('assemblies'))) :
                 numpy_table[i,a] = chunk.count(str(a))
        # print lines
        for a in range(len(config.get('assemblies'))) :
            line = f"{loc}\t{name}\t"+str(a+1)+"\t"
            for i in range(0, int(len(seq)/slice_size)+1):
                line = line+str(numpy_table[i,a])+","
            line = line[:-1]
            #print(line)
            slices.append(line)

    #print(lines)  
    return slices
    
# Function to process assemblies and store sequence names and lengths
def process_assemblies(config):
    assemblies = config.get('assemblies', [])
    all_sequences = {}

    for assembly_file in assemblies:
        sequences = parse_fasta(assembly_file)
        all_sequences[assembly_file] = sequences
        #print(f"Processed assembly '{assembly_file}':")
        #for seq_name, seq_len in sequences:
        #    print(f"  Sequence: {seq_name}, Length: {seq_len} bp")

    return all_sequences

# Function to validate the intervals
def check_intervals(config):
    intervals = config.get('intervals', [])
    previous_end = None  # To keep track of the end value of the previous interval

    for i, interval in enumerate(intervals):
        try:
            # Split the interval by comma and convert to floats
            start, end = map(float, interval.split(','))
            
            # Check that the second value is larger than the first
            if start >= end:
                print(f"Invalid interval: '{interval}' (the second value must be larger than the first).", file=sys.stderr)
                sys.exit(1)  # Stop execution on failure
            else:
                print(f"Interval '{interval}' is valid.")
            
            # Check that the last value of the previous interval is smaller than the current start
            #print(previous_end, start)
            if previous_end is not None and previous_end > start:
                print(f"Invalid order: '{interval}' (the start value must be larger than the end of the previous interval).", file=sys.stderr)
                sys.exit(1)  # Stop execution on failure

            # Update the previous_end value for the next comparison
            previous_end = end

        except ValueError:
            # Catch any non-numeric values and report an error
            print(f"Invalid interval: '{interval}' (must be numeric values).", file=sys.stderr)
            sys.exit(1)  # Stop execution on failure

# Function to check that the number of intervals equals the number of assemblies
def check_intervals_and_assemblies_count(config):
    intervals = config.get('intervals', [])
    assemblies = config.get('assemblies', [])

    if len(intervals) != len(assemblies):
        print(f"Mismatch: Number of intervals ({len(intervals)}) does not match number of assemblies ({len(assemblies)}).", file=sys.stderr)
        sys.exit(1)  # Stop execution on failure
    else:
        print(f"Number of intervals matches the number of assemblies: {len(intervals)}")

# Function to check if all specified software can be run on the system
def check_softwares(config):
    softwares = config.get('softwares', [])
    for software in softwares:
        try:
            # Run the command with no arguments
            result = subprocess.run([software], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            # print(result.returncode)       
            # Check the return code
            if result.returncode == 0 or result.returncode == 1 :
                print("The "+software+" command is accessible.")
            else:
                print("The "+software+" command is not accessible.")
                sys.exit(1) 
                return False
        except FileNotFoundError:
            print("The "+software+" command is not found in the environment.")
            return False
        except Exception as e:
            print(f"An error occurred: {e}")
            return False

def create_kmc_tmp_directory():
    #print("create_kmc_tmp_directory")
    try:
        result = subprocess.run("mkdir tmp", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode == 0 :
            print("tmp directory created")
        else:
            print("already existing tmp directory")
    except FileNotFoundError:
        print("the tmp directory could not be created")
        return False

def generate_kmc_bases_per_interval(config):
    intervals = config.get('intervals', [])
    reads_f_path = config.get('reads')
    reads_format = config.get('kmc_reads_file_format')
    kmc_memory = config.get('kmc_memory')
    kmc_cpus = config.get('kmc_cpus')
    prefix = config.get('prefix')
    kmer_length = config.get('kmer_length')
    
    for i, interval in enumerate(intervals):
        start, end = map(str,interval.split(','))
        kmc_line  =  f"kmc -v -k{kmer_length} -m{kmc_memory} -ci{start} -cx{end} -cs10000 -t{kmc_cpus} -{reads_format} @{reads_f_path} {prefix}.{start}-{end}.k{kmer_length} tmp"
        print(kmc_line)
        try:
            result = subprocess.run([kmc_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"{prefix}.{start}-{end}.k{kmer_length} database created")
            else:
                print(f"{prefix}.{start}-{end}.k{kmer_length} database not created")
        except FileNotFoundError:
            print(f"{prefix}.{start}-{end}.k{kmer_length} database could not be created")
            return False

def dump_kmc_bases_per_interval(config):
    intervals = config.get('intervals', [])
    reads_f_path = config.get('reads')
    reads_format = config.get('kmc_reads_file_format')
    kmc_memory = config.get('kmc_memory')
    kmc_cpus = config.get('kmc_cpus')
    prefix = config.get('prefix')
    kmer_length = config.get('kmer_length')    
    # max_interval = intervals[len(intervals)-1]
    #print("intervals ", intervals[1], len(intervals), max_interval)

    # dumping kme database in file
    for i, interval in enumerate(intervals):
        start, end = map(str,interval.split(','))
        kmc_line  =  f"kmc_dump {prefix}.{start}-{end}.k{kmer_length} {prefix}.{start}-{end}.k{kmer_length}.dump"
        print(kmc_line)
        try:
            result = subprocess.run([kmc_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"{prefix}.{start}-{end}.k{kmer_length} database created")
            else:
                print(f"{prefix}.{start}-{end}.k{kmer_length} database not created")
        except FileNotFoundError:
            print(f"{prefix}.{start}-{end}.k{kmer_length} database could not be created")
            return False

    # transforming dump file in fasta with read count = interval rank   
    fo = open(f"all.k{kmer_length}.dump.fa","w")
    for i, interval in enumerate(intervals):
        start, end = map(str,interval.split(','))
        fi = open(f"{prefix}.{start}-{end}.k{kmer_length}.dump","r")
        for ind, line in enumerate(fi) :
            kmer, count = map(str,line.split('\t'))
            for ii in range(i+1) :
                fo.write(">"+str(ind)+"-"+str(ii)+"\n")
                fo.write(kmer+"\n")
        fi.close()
    fo.close()
    last_dump_file = f"{prefix}.{start}-{end}.k{kmer_length}.dump"
    
    # find kmers in assemblies 
    # create kmer jellyfish database 
    jellyfish_line  =  f"jellyfish count -m {kmer_length} -C -o All.jf -s 100M -t 10 all.k{kmer_length}.dump.fa"
    print(jellyfish_line)
    try:
        result = subprocess.run([jellyfish_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode == 0 :
            print("Jellyfish database created")
        else:
            print("Jellyfish database not created")
    except FileNotFoundError:
        print("Jellyfish database could not be created")
        return False
    
    # searching all the assemblies with this data base 
    locate = []
    assemblies = config.get('assemblies', [])
    for i, assembly in enumerate(assemblies):
        jellyfish_line  =  f"query_per_sequence All.jf {assembly} > {assembly}.locate"
        print(jellyfish_line)
        locate.append(f"{assembly}.locate")
        try:
            result = subprocess.run([jellyfish_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"Jellyfish kmer locate for {assembly}")
            else:
                print(f"Jellyfish kmer not locate for {assembly}")
        except FileNotFoundError:
            print(f"Jellyfish location did not work for {assembly}")
            return False        
    return locate, last_dump_file     

# produces all the links with a block
def process_kmer_link_block(block) :
    res =  []
    b = []
    if len(block) > 1 : 
        for i, l in block.items() :
            contig, pos, kmer = map(str,l.split('\t'))
            b.append(contig)
        sb = list(set(b))
        #print(sb)
        if len(sb) > 1 :
            res = [f"{a},{b}" for idx, a in enumerate(sb) for b in sb[idx + 1:]]
            #print("res", res)
    return res

def split_wrong_line(line, kmer_length):
    lines = []
    b = line.split("\t")
    #print("wrong line :", line)
    # correct split contatenated lines
    for i in range(len(b) -3, 0, -2) :
        #print(i, "bi : " ,b[i])
        c=b[i][0:kmer_length]
        d = b[i][kmer_length:]
        b = b[:i] + [c] + [d] + b[i + 1:]
        #print(b)

    # rebuild lines
    for i in range(0,len(b),3) :
        lines.append(b[i]+"\t"+b[i+1]+"\t"+b[i+2])

    #print("corrected lines", lines)

    return lines

    
def locate_common_kmers_in_contig_pairs(config, last_dump_file):
    intervals = config.get('intervals', [])
    prefix = config.get('prefix')
    kmer_length = config.get('kmer_length')  
        
    # index all assemblies with samtools faidx
    assemblies = config.get('assemblies', [])
    for i, assembly in enumerate(assemblies):
        samtools_line  =  f"samtools faidx {assembly}"
        try:
            result = subprocess.run([samtools_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"samtools faidx {assembly}")
            else:
                print(f"samtools faidx failed {assembly}")
        except FileNotFoundError:
            print(f"samtools faidx did not work for {assembly}")
            return False            

    #get kmer file for the top interval 
    fi = open(last_dump_file,"r")
    fo = open(f"{last_dump_file}.first_column","w")
    for ind, line in enumerate(fi) :
        kmer, count = map(str,line.split('\t'))
        fo.write(kmer+"\n")
    fi.close()   
    fo.close() 
    
    #and locate these kmer in each assembly with seqc 
    for i, assembly in enumerate(assemblies):
        seq_line  =  f"seqc run find_kmer.seq {assembly} {last_dump_file}.first_column > {prefix}.{assembly}.localisation"
        print(seq_line)
        try:
            result = subprocess.run([seq_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"search kmers with seqc {assembly}")
            else:
                print(f"search kmers with seqc failed {assembly}")
        except FileNotFoundError:
            print(f"search kmers with seqc did not work for {assembly}")
            return False            

    # cleaning localisation for problematic lines (empty or duplicates)
    assemblies = config.get('assemblies', [])
    for i, assembly in enumerate(assemblies):
        fi = open(f"{prefix}.{assembly}.localisation","r")
        fo = open(f"{prefix}.{assembly}.localisation.clean","w")
        block = []
        ll = 0
        for ind, line in enumerate(fi) :
            if len(line.split('\t')) == 3 :
                fo.write(line)
            else :
                if len(line.split('\t')) > 3 :
                    # split line with function and return an array of correct lines 
                    lino = split_wrong_line(line[:-1], kmer_length)
                    for li in lino :
                        fo.write(li+"\n")

        fi.close()   
        fo.close()     
                
    # sort kmer localisation by profile 
    assemblies = config.get('assemblies', [])
    for i, assembly in enumerate(assemblies):
        sort_line  =  f"sort -k3,3 {prefix}.{assembly}.localisation.clean > {prefix}.{assembly}.localisation.clean.sort "
        try:
            result = subprocess.run([sort_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"sort kmer localisation {assembly}")
            else:
                print(f"sort kmer localisation {assembly}")
        except FileNotFoundError:
            print(f"sort kmer localisation did not work for {assembly}")
            return False                 

    # create kmer couples 
    for i, assembly in enumerate(assemblies):
        fi = open(f"{prefix}.{assembly}.localisation.clean.sort","r")
        fo = open(f"{prefix}.{assembly}.localisation.clean.sort.links","w")
        block = {}
        for ind, line in enumerate(fi) :
            #print(line)
            contig, pos, kmer = map(str,line.split('\t'))
            if ind != 0 :
                if kmer == refkmer :
                    block[contig]= line
                else :
                    if len(block) > 1 :
                        # getting all the links between contigs
                        res = process_kmer_link_block(block)
                        # here we should print the couples of kmer information and not only the link between the contigs
                        for r in res:
                            r1 = r.split(",")
                            #print(block[r1[0]][:-1]+"\t"+block[r1[1]][:-1])
                            if block[r1[0]][2] > block[r1[1]][2] :
                                fo.write(block[r1[1]][:-1]+"\t"+block[r1[0]][:-1]+"\n")
                            else :
                                fo.write(block[r1[0]][:-1]+"\t"+block[r1[1]][:-1]+"\n")
                    block = {}
                    block[contig]= line
                    refkmer = kmer
            else :
                block[contig]= line
                refkmer = kmer
            
        fi.close()   
        fo.close()     

def select_contig_pairs(config):
    prefix = config.get('prefix')
    assemblies = config.get('assemblies', [])
    
    # for all haplotypes find contig pair list with links counts
    for i, assembly in enumerate(assemblies):
        count_line  =  f"cut -f1,4 {prefix}.{assembly}.localisation.clean.sort.links | sed 's/\t/,/' | sort -k1,1 | uniq -c > {prefix}.{assembly}.localisation.clean.sort.links.list"
        try:
            result = subprocess.run([count_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if result.returncode == 0 :
                print(f"count kmer couples {assembly}")
            else:
                print(f"count kmer couples failed {assembly}")
        except FileNotFoundError:
            print(f"count kmer couples did not work for {assembly}")
            return False

# for each haplotype generate de profile for both contigs of the pair : for the complete kmers and the corresponding kmers 
def generate_count_files_for_graphs(config, count_lines):
    prefix = config.get('prefix')
    assemblies = config.get('assemblies', [])
    slice_size = config.get('slice_size')
    minimum_link_count = config.get('minimum_link_count')
    contig_chunks = {}
    
    # read pair file 
    for i, assembly in enumerate(assemblies):
        contigs = []
        contig_pairs = []
        fi = open(f"{prefix}.{assembly}.localisation.clean.sort.links.list","r")
        fi2 = open(f"{prefix}.{assembly}.localisation.clean.sort.links","r")
        fo = open(f"{prefix}.{assembly}.localisation.clean.sort.links.list.general_stats","w")
        fo2 = open(f"{prefix}.{assembly}.localisation.clean.sort.links.list.pair_stats","w")
        for ind, line in enumerate(fi) :
             line = '\t'.join([x for x in line[:-1].split(' ') if len(x)>0])
             contigs.append(line.split("\t")[1].split(",")[0])
             contigs.append(line.split("\t")[1].split(",")[1])
             contig_pairs.append([line.split("\t")[1].split(",")[0],line.split("\t")[1].split(",")[1]])
        #print("contigs for which to get general stats", set(contigs))
        #print("contig pair for wich to get local stats", contig_pairs)
        
        for i, line in enumerate(count_lines) :
            #print(line)
            if line.split("\t")[0] == assembly :
                if line.split("\t")[1] in contigs :
                    # add contig chunks to dictionary
                    if not line.split("\t")[1] in contig_chunks : 
                        contig_chunks[line.split("\t")[1]] = len(line.split("\t")[3].split(","))
                    # write stat line
                    fo.write(line+"\n")
        
        #print("contig chunks", contig_chunks)
        counts = {}          
        for ind, line in enumerate(fi2) :
            #print(line)
            for pair in contig_pairs :
                # if link matches to contig pair  
                if (line.split("\t")[0] == pair[0] and line.split("\t")[3] == pair[1]) or (line.split("\t")[3] == pair[0] and line.split("\t")[0] == pair[1]) :
                   # put the pair name in the same order for both pairs 
                   if line.split("\t")[0] > line.split("\t")[0] :
                       name = line.split("\t")[0]+"-"+line.split("\t")[3]
                   else :
                       name = line.split("\t")[3]+"-"+line.split("\t")[0]
                   if name+"\t"+line.split("\t")[0]+"\t"+str(int(int(line.split("\t")[1])/slice_size)) in counts :
                       counts[name+"\t"+line.split("\t")[0]+"\t"+str(int(int(line.split("\t")[1])/slice_size))] += 1
                   else :
                       counts[name+"\t"+line.split("\t")[0]+"\t"+str(int(int(line.split("\t")[1])/slice_size))] = 1
                   if name+"\t"+line.split("\t")[3]+"\t"+str(int(int(line.split("\t")[4])/slice_size)) in counts :
                       counts[name+"\t"+line.split("\t")[3]+"\t"+str(int(int(line.split("\t")[4])/slice_size))] += 1
                   else :
                       counts[name+"\t"+line.split("\t")[3]+"\t"+str(int(int(line.split("\t")[4])/slice_size))] = 1
                       
        # Dictionary to store the sorted counts per group/name combination
        grouped_data = defaultdict(lambda: defaultdict(dict))

        # Step 1: Parse the keys and organize by group, name, and position
        for key, count in counts.items(): 
            group, name, position = key.split('\t')
            position = int(position)  # Convert position to integer for sorting
            grouped_data[group][name][position] = count

        # print(grouped_data)
        # Step 2: Sort by position and prepare the final array
        result = []

        for group, names in grouped_data.items():
            for name, positions in names.items():
                # Sort the counts by position
                sorted_positions = sorted(positions.keys())
                sorted_counts = [count for pos, count in sorted(positions.items())]
                result.append([group, name, sorted_positions, sorted_counts])

        # Print the final array
        #print(result)
        #print(contig_chunks)
        #print(result[1], result[1][1])
        for i in result :
            #print(i)
            #print(contig_chunks[i[1]])
            numpy_table = np.zeros(contig_chunks[i[1]])
            #print(numpy_table)
            for k, j in enumerate(i[2]) :
                numpy_table[j] = i[3][k]
            #print(numpy_table)
            line = assembly+"\t"+i[0]+"\t"+i[1]+"\t"
            for i in range(contig_chunks[i[1]]):
                line = line+str(numpy_table[i])+","
            line = line[:-1]
            fo2.write(line+"\n")
        
        fi.close()
        fo.close()
        fo2.close()
        
# Main function to handle command-line arguments and run the checks
def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Running manual purging pipeline to find contigs with duplicated kmers expected only once in the assembly")
    
    # Add argument for the configuration file path
    parser.add_argument('config_file', type=str, help="Path to the YAML configuration file")
    
    # Parse command-line arguments
    args = parser.parse_args()

    # Load the configuration from the specified file
    config = load_yaml_config(args.config_file)

    # Check if the reads and assemblies files are accessible
    check_files(config)

    # Validate the intervals
    check_intervals(config)

    # Check if the number of intervals equals the number of assemblies
    check_intervals_and_assemblies_count(config)

    # Check if all the specified software can be run
    check_softwares(config)

    # create temporary directory needed for kmc
    create_kmc_tmp_directory()
    
    # Process the assemblies and store sequence names and lengths
    if int(config.get('starting_step')) <= 1 :
        process_assemblies(config)

    # produce kmc interval databases
    if int(config.get('starting_step')) <= 2 :
        generate_kmc_bases_per_interval(config)
    
    # produce kmc interval databases
    if int(config.get('starting_step')) <= 3 :
        locate, last_dump_file = dump_kmc_bases_per_interval(config)  
        f = open('file1.pic', 'wb')  
        pickle.dump([locate, last_dump_file], f)
        f.close()
    
    # process query_per_sequence count file / BEWARE this should be produced only for selected contigs not for all 
    if int(config.get('starting_step')) <= 4 :
        try:
            locate
        except NameError:
            f = open('file1.pic', 'rb')
            [locate, last_dump_file] = pickle.load(f)
            print("retrieved from pickle file1 : ",locate)
            f.close()
        slice_count = []
        for i, loc in enumerate(locate) :
            counts = parse_count(loc)
            slice_count = count_per_slice(loc, counts,config, slice_count)
        f = open('file2.pic', 'wb')  
        pickle.dump([counts, slice_count], f)
        f.close()
    #print(slice_count)
    
    # find duplicated kmers in all assemblies using the last kmer : 
    if int(config.get('starting_step')) <= 5 :
        try:
            last_dump_file
        except NameError:
            f = open('file1.pic', 'rb')
            [locate, last_dump_file] = pickle.load(f)
            print("retrieved from pickle file1 : ",locate)
            f.close()
        locate_common_kmers_in_contig_pairs(config, last_dump_file)
    
    # find contig links to select contig couples to manually purge 
    if int(config.get('starting_step')) <= 6 :
        select_contig_pairs(config)

    # get kmer location for each pair of contigs / the links are in localisation.clean.sort.links
    if int(config.get('starting_step')) <= 7 :
        try:
            slice_count
        except NameError:
            f = open('file2.pic', 'rb')
            [counts, slice_count] = pickle.load(f)
            #print("retrieved from pickle file2 : ",slice_count)
            f.close()    
        generate_count_files_for_graphs(config, slice_count)
    
if __name__ == "__main__":
    main()

