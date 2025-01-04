import pandas as pd
import re

def convert_to_tuples(input_string):
    # Use a regular expression to find all occurrences of a number followed by a letter
    pattern = r'(\d+)([A-Z])'
    
    # Find all matches and convert them into tuples (number, symbol)
    result = [(int(num), symbol) for num, symbol in re.findall(pattern, input_string)]
    # print("List of tuple ",result)
    return result

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

bam_path = './aligned.bam'  # Replace with your BAM file path
df = pd.read_csv(bam_path, sep='\t',skiprows=2, header=None, on_bad_lines='skip')
# df = bam_to_dataframe(bam_path)
df =df[[0,3,5,9,10]]
df.rename(columns={0:"Id",3:"Start_Index",5:"CIGAR",9:"Read",10:"Fred-Score"}, inplace=True)
# print(df.head())/
df_copy = df
for ind,row in df_copy.iterrows():
    new_str, prcessed_fred = process_CIGAR(row)
    
    # break
    df.loc[ind,"Read"] = new_str
    df.loc[ind,"Fred-Score"] = ','.join(str(x) for x in prcessed_fred)
    # print(len(new_str))
    # print(len(prcessed_fred))
    # print(len(prcessed_fred))
    # break
df = df[["Start_Index","Read","Fred-Score"]]
# df = df[["Id","Start_Index","Read","Fred-Score"]]
# df = df[["Start_Index","Read"]]
df["Start_Index"] = df["Start_Index"]-1
df.to_csv("processed_bam_data.csv",index = False, sep = '|')