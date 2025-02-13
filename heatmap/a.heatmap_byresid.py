import os
import sys
from dotenv import dotenv_values

if len(sys.argv) != 4:
    print("Usage: python3 script.py alignerInputs.txt simReplicates.txt nagResidues.txt")
    sys.exit(1)

alignerInputs = sys.argv[1]
simReplicates = sys.argv[2]
nagResidues = sys.argv[3]

config = dotenv_values(alignerInputs)
globals().update(config)

replicates = []
with open(simReplicates, "r") as simReplicates:
    for line in simReplicates:
        line = line.strip()
        replicates.append(line)

nag_residues = []
with open(nagResidues, "r") as nagResidues:
    for line in nagResidues:
        line = line.strip()
        nag_residues.append(line)
print(nag_residues)

def extract_columns(filename):
    clash_list = [] 
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(PREFIX):     # Lines in chimera output file that start with prefix contain overlap info
                columns = line.strip().split()
                system_residue_in_clash = columns[-4]   # Column in chimera file with residue number of clash
                clash_list.append(system_residue_in_clash)
    return clash_list

for rep in replicates:
    heatmap_directory_by_rep = f"{rep}-heatmap"
    if os.path.exists(heatmap_directory_by_rep):
        pass
    else:
        os.mkdir(heatmap_directory_by_rep)
        print(f"{heatmap_directory_by_rep} created ...")
    TEMPLATE_PATH_TO_OVERLAPS = TEMPLATE_PATH_TO_OVERLAPS.format(REP=f"{rep}-overlaps")
    print(TEMPLATE_PATH_TO_OVERLAPS)
    residue_list = []
    clash_record = []
    system_filename = SYSTEM
    with open(system_filename, 'r') as system_file:
        system_lines = system_file.readlines()
    for line in system_lines:
        if line.startswith('ATOM'):
            elements = line.split()
            # excluding water, Na+ atoms
            if not 'WAT'in elements and not 'Na+' in elements:
                resid_num = elements[4]
                residue_list.append(resid_num)
                clash_record.append([resid_num])
    
    print('Number of atoms:',len(residue_list))
    # print(residue_list)
    for nag in nag_residues:
        heatmap_directory_by_nag = f"{nag}.heatmap"
        heatmap_final_path_by_rep_and_nag = os.path.join(heatmap_directory_by_rep,heatmap_directory_by_nag)
        if os.path.exists(heatmap_final_path_by_rep_and_nag):
            pass
        else:
            os.mkdir(heatmap_final_path_by_rep_and_nag)
            print(f"{heatmap_final_path_by_rep_and_nag} created ...")
        nag_directory = f"{nag}.outputs"
        fileCount = 1
        while fileCount <= int(STRUCTURE_LIMIT_CHIMERA):
            overlap_file = f"{rep}-ss{fileCount}-proteinOverlaps.txt"
            overlap_path = os.path.join(TEMPLATE_PATH_TO_OVERLAPS,nag_directory,overlap_file)
            clash = extract_columns(overlap_path)
            fileCount += 1
        file_record = [] 
        for resid in residue_list:
            if resid in clash:
                file_record.append(1)
            else:
                file_record.append(0)

        for i in range(len(residue_list)):
            clash_record[i].append(file_record[i])
        
        
        B_val = [(sum(c_record[1:])/int(STRUCTURE_LIMIT_CHIMERA)) * 100 for c_record in clash_record]
        print('length of B_val',len(B_val))     #The length of B_val should match the length of atom_list
        #write B_val to a file:
        print(f"Writing B values to {heatmap_final_path_by_rep_and_nag}/Bval_byresid.csv...")
        with open(f'{heatmap_final_path_by_rep_and_nag}/Bval_byresid.csv', 'w') as file:
            for val in B_val: 
                file.write(str(val) + '\n')
        print("...completed")
       
        # Generate an updated version of SYSTEM.pdb that contains updated B-factors
        
        modified_lines = []
        with open(system_filename, "r") as system_file:
            for line in system_file:
                elements = line.split()
                cond1 = line.startswith("ATOM")
                if cond1:
                    cond2 = 'WAT' not in elements
                    cond3 = 'Na+' not in elements
                if cond1 and cond2  and cond3:
                    b_value = B_val.pop(0)
                    modified_line = line[:60] + "{:>6.2f}".format(b_value) + line[66:]
                    modified_lines.append(modified_line)
                else:
                    modified_lines.append(line)
        modified_pdb_filename = f"modified_B-factor_byresid.pdb"
        with open(f"{heatmap_final_path_by_rep_and_nag}/{modified_pdb_filename}", 'w') as final:
            final.writelines(modified_lines)
        print(f"...Updated pdb file: {heatmap_final_path_by_rep_and_nag}/{modified_pdb_filename}")