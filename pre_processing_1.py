import pandas as pd
import re

def convert_to_tuples(input_string):
    # Use a regular expression to find all occurrences of a number followed by a letter
    pattern = r'(\d+)([A-Za-z])'
    
    # Find all matches and convert them into tuples (number, symbol)
    result = [(int(num), symbol) for num, symbol in re.findall(pattern, input_string)]
    # print("List of tuple ",result)
    return result

def process_CIGAR(row):
    list_of_tuples = convert_to_tuples(row['CIGAR'])
    curr_str = row['Read']
    # l = len(curr_str)
    # print(f'Current String:len{l} ',curr_str)
    index = -1
    new_str = ''
    for x,symbol in list_of_tuples:
        if symbol == 'M':
            new_str += curr_str[index+1:index+x+1]
        if symbol == 'I':
            pass
        if symbol == 'D':
            new_str += '_'*x
        index += x 
    # l =len(new_str)
    # print(f"Processed String {l}: ",new_str)
    return new_str


bam_path = '/home/tarunpal/Desktop/temp/Mayank_Project/aligned.bam'  # Replace with your BAM file path
df = pd.read_csv(bam_path, sep='\t',skiprows=2, header=None, on_bad_lines='skip')
# df = bam_to_dataframe(bam_path)
df =df[[0,3,5,9,10]]
df.rename(columns={0:"Id",3:"Start_Index",5:"CIGAR",9:"Read",10:"Fred-Score"}, inplace=True)
# print(df.head())/
df_copy = df
for ind,row in df_copy.iterrows():
    new_str= process_CIGAR(row)
    df.loc[ind,"Read"] = new_str
    # print(new_str)
    # break
df = df[["Id","Start_Index","Read","Fred-Score"]]
df.to_csv("processed_bam_data.csv",index = False)