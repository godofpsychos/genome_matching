import pandas as pd
import re

def convert_to_tuples(input_string):
    # Use a regular expression to find all occurrences of a number followed by a letter
    pattern = r'(\d+)([A-Z])'
    
    # Find all matches and convert them into tuples (number, symbol)
    result = [(int(num), symbol) for num, symbol in re.findall(pattern, input_string)]
    # print("List of tuple ",result)
    return result
def get_second_line(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if len(lines) >= 2:
                return lines[1].strip()  # Return the second line (1-based index)
            else:
                return ""
    except FileNotFoundError:
        print("Could not open the file!")
        return ""
def process_CIGAR(row):
    list_of_tuples = convert_to_tuples(row['CIGAR'])
    curr_str = row['Read']
    prcessed_fred = get_quality_score(row["Fred-Score"])
    l = len(curr_str)
    # print(l)
    # print(f'Current String:len{l} ',curr_str)
    index = -1
    fred = []
    new_str = ''
    for x,symbol in list_of_tuples:
        if symbol == 'M':
            new_str += curr_str[index+1:index+x+1]
            fred.extend(prcessed_fred[index+1:index+x+1])
        if symbol == 'I':
            pass
        if symbol == 'D':
            new_str += '_'*x
            fred.extend([-1]*x)
        index += x 
    # l =len(new_str)
    # print(f"Processed String {l}: ",new_str)
    return new_str,fred

def get_quality_score(row):
    lis = [10 ** (-1 *((ord(x)-33)/10 ))for x in row]
    # print(len(lis))
    # lis = ','.join(str(x) for x in lis)
    return lis



example={
    "index": 0,
    "Unmatch": 0,
    "Match": 0,
    "Fred-prob":[]
}
prcessed_fred =[]
bam_path = './aligned.bam'  # Replace with your BAM file path
ref_gnome_file_path = "/home/tarun/Desktop/genome_matching/EPI_ISL_402124-ORF1ab.fasta"
ref_gnome = get_second_line(ref_gnome_file_path)
    
df = pd.read_csv(bam_path, sep='\t',skiprows=2, header=None, on_bad_lines='skip')
# df = bam_to_dataframe(bam_path)
df =df[[0,3,5,9,10]]
df.rename(columns={0:"Id",3:"Start_Index",5:"CIGAR",9:"Read",10:"Fred-Score"}, inplace=True)
# print(df.head())/
df_copy = df
col_dict = {}
for ind,row in df_copy.iterrows():
    new_str, prcessed_fred = process_CIGAR(row)
    
    # break
    df.loc[ind,"Read"] = new_str
    df.loc[ind,"Fred-Score"] = ','.join(str(x) for x in prcessed_fred)
    start = row["Start_Index"]
    read = new_str    
    fred = prcessed_fred    
    if start < 0:
        continue
    for i in range(start, min(len(ref_gnome),start+len(read))):
        if i not in col_dict:
            col_dict[i] = {
                            "Unmatch": 0,
                            "Match": 0,
                            "Fred-prob":[]}
        # print(i,[i-start],ref_gnome[i])
        col_dict[i]["Match"] += 1 if read[i-start] == ref_gnome[i] else 0
        col_dict[i]["Unmatch"] += 1 if read[i-start] != ref_gnome[i] else 0
        col_dict[i]["Fred-prob"].append(fred[i-start])

for idx, value in col_dict.items():
    value["Fred-prob"] = ','.join(map(str, value["Fred-prob"]))  # Convert list to string

# Convert col_dict to a pandas DataFrame
new_df = pd.DataFrame.from_dict(col_dict, orient='index')
df = df[["Start_Index","Read","Fred-Score"]]
# df = df[["Id","Start_Index","Read","Fred-Score"]]
# df = df[["Start_Index","Read"]]
df["Start_Index"] = df["Start_Index"]-1
# for idx,row in df.iterrows():
#     start = row["Start_Index"]
#     read = row["Read"]
#     fred = row["Fred-Score"]

df.to_csv("processed_bam_data.csv",index = False, sep = '|')
new_df.to_csv("bam_data_summary.csv",index = True, sep = '|')