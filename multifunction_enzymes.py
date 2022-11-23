__usage__="""
                    multifunction_enzymes.py
                    --in <INPUT DIRECTORY: KIPES result directory>
                    --out <OUTPUT_FILE>
                    --enzymes <ENZYMES_OF_INTEREST>
                  
                    """

import sys,glob, pandas as pd
# --- end of imports --- #
   
def ext_cons_res (filename1,filename2):
    """Extracts enzyme candidate genes with perfectly conserved residues from the conserved_residues.txt files"""
    pref_seq=[]
    with open(filename1, "r") as f1:
        line=f1.readline()
        while line:
            parts=line.strip().split('\t')
            length=len(parts)
            if length >=2:  #executes the code only if there is residues data in the file
                ele=parts[1]
                chk= True
                for item in parts[1:]:  #finds the candidate enzymes with perfect conserved residues 
                    if ele!=item:
                        chk= False
                    if parts[1]=="-" or parts[-1]=="-":
                        if "False" in line:
                            chk=False
                            break;
                if (chk==True):
                    pref_seq.append(parts[0])
                line=f1.readline()
    """"Extract the corresponding peptide sequences for the enzyme candidates from the final_pep_files(.fasta)"""
    result = []
    lines = open(filename2, 'r').read().splitlines()
    for i in range(len(lines)- 1): #read all the lines from the line
        if lines[i][1:] in pref_seq: #if header matches the enzyme candidate
            result.append((lines[i][1:], lines[i+1])) #returns the enzyme ID and its sequence
    return result
         
def main(arguments):
    input_folder= arguments[ arguments.index('--in')+1 ]
    output_file=arguments[arguments.index('--out')+1]
    enzymes_of_interest=arguments[arguments.index('--enzymes')+1]
    if input_folder[-1] != "/":
        input_folder+="/" 
    if "," in enzymes_of_interest:
        enzymes_of_interest=enzymes_of_interest.split(',')
    else:
        enzymes_of_interest=[ enzymes_of_interest ]
    cols = ["Species", "Enzyme", "ID", "Sequences"]
    data = []
    for enzyme in enzymes_of_interest:
        cons_res_files=glob.glob(input_folder+ "*/conserved_residues/"+enzyme+"_conserved_residues.txt", recursive=True)
        seq_files=glob.glob(input_folder+ "*/final_pep_files/"+enzyme+".fasta", recursive=True)
        for filename1, filename2 in zip(cons_res_files, seq_files): #takes input file in the same order
            cands_seqs = ext_cons_res(filename1,filename2)
            spec=filename1.split('/')[-3]
            ide=filename1.split('/')[-1]
            ids=ide.split('_')[0]
            for cand, seq in cands_seqs:
                    data.append([spec, ids, cand, seq])
    df = pd.DataFrame(data, columns=cols)
    df = df.sort_values('Species')  #sort columns by species name
    IDs = df.ID.unique()
    groups = [] 
    for i in IDs:
        groups.append(df.groupby('ID').get_group(i))  #Groups the Dataframe on the basis of IDs
    final_groups = []
    for i in range(len(groups)):
        if groups[i].shape[0] > 1: 
            final_groups.append(groups[i])  #only show those rows where more than 1 enzyme candidates match
    if len(final_groups) != 0:   #no enzyme candidates match
        final_df = pd.concat(final_groups).to_csv(output_file, sep='\t', index=False )
        
    else: 
        print('No common enzyme candidates found. Cannot prove enzyme multifunctionality.')
                   
            
if '--in' in sys.argv and '--out' in sys.argv and '--enzymes' in sys.argv:
    main( sys.argv )
else:
    sys.exit(__usage__)

                        