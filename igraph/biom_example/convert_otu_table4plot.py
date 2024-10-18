import pandas as pd
import re

def taxonomy_parser(taxonomy_str,db_type):
    
    if db_type == "Silva":
        separator = ';'
    elif db_type == "GreenGenes":
        separator = r"\|"
    else:
        print(db_type,"不是Silva或者GreenGenes")
        return False
    '''

    # kingdom,phylum,class,order,family,genus,species
    #  0     , 1   , 2   ,3     ,  4     ,5   ,6
    #  1        2      3     4     5      6     7

    D_0__Bacteria; D_1__Firmicutes; D_2__Clostridia; D_3__Clostridiales; D_4__Christensenellaceae; D_5__Christensenellaceae R-7 group
    D_0__Bacteria; D_1__Bacteroidetes; D_2__Bacteroidia; D_3__Bacteroidales; D_4__Muribaculaceae; D_5__uncultured bacterium; D_6__uncultured bacterium
    '''
    records = re.split(separator,taxonomy_str.strip()) # 未知长度,可能不完整
    taxonomy_type_idx_dic = {
        'kingdom':0,
        'phylum':1,
        'class':2,
        'order':3,
        'family':4,
        'genus':5,
        'species':6
    }
    taxonomy_lst = [] #结果
    for index in range(len(taxonomy_type_idx_dic)):
        if index < len(records):
            # 根据数据库内容切换，正则匹配切换
            if db_type == "Silva":
                pattern = re.compile(r"D_"+str(index)+"__(.+)")
            elif db_type == "GreenGenes":
                pattern = re.compile(r"[kpcofgs]__(.+)")
            match_obj = re.search(pattern,records[index].strip())
            if match_obj :
                taxonomy_info = match_obj.group(1)
                # 特殊符号替换下划线_
                if "[" in taxonomy_info:
                    taxonomy_info = re.split("[\[\]]",taxonomy_info)[1]
                taxonomy_info = re.sub('\s','_',taxonomy_info.strip())
                # 结果修改
                #taxonomy_info = "D_"+str(index)+"__"+taxonomy_info
                taxonomy_lst.append(taxonomy_info)
            else: # 未能匹配
                #print("Warring!：not match:\t\""+records[index]+"\"\n")
                taxonomy_lst.append('NA')
        else: # 超出范围
            taxonomy_lst.append('NA')
    return taxonomy_lst



# taxonomy_str="D_0__Bacteria; D_1__Firmicutes; D_2__Clostridia; D_3__Clostridiales; D_4__Christensenellaceae; D_5__Christensenellaceae R-7 group"
# taxonomy_str2=" k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_anginosus_group"

# print(taxonomy_parser(taxonomy_str,db_type="Silva"))
# print(taxonomy_parser(taxonomy_str2,db_type="GreenGenes"))


df = pd.read_csv('otu_table.txt', sep='\t', skiprows=1,encoding='utf-8')


out_table_df= df.drop('taxonomy', axis=1, inplace=False)
out_table_df = out_table_df.T

# 转置过，所以 header = False ；使用第一行作为列
out_table_df.to_csv('otu_table.tsv',sep='\t',index=True,header=False)

# 创建新的列来存储解析后的分类信息
new_columns = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for col in new_columns:
    df[col] = None

# 遍历 'taxonomy' 列并使用 taxonomy_parser 函数解析
for index, row in df.iterrows():
    if pd.notna(row['taxonomy']):
        parsed_taxonomy = taxonomy_parser(row['taxonomy'], db_type="Silva")  # 根据实际情况选择数据库类型
        for i, col in enumerate(new_columns):
            # print(i, col) 迭代每个 新列名，将每行对应的新列明 单个元素进行修改
            df.at[index, col] = parsed_taxonomy[i]

# 删除原始的 'taxonomy' 列
df.drop('taxonomy', axis=1, inplace=True)
df.to_excel("taxonomy_result.xlsx",index=False)

########################## 只保留分类表 ############################
index_sample = len(df.columns.to_list()) -7
selected_columns_4_drop = df.columns[1:index_sample].to_list()
df.drop(selected_columns_4_drop, axis=1, inplace=True)
df.to_csv('otu_pro.tsv',sep='\t',index=False)
##################################################################
