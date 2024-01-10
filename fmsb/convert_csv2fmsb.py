import  re
import pandas as pd


MAX = 0.5
MIN = 0.0

# 读取目标列信息
target_genus_dic = {}
with open('20_Genus.txt',mode='rt',encoding='utf-8') as fh:
    for line in fh:
        key_raw = line.strip()
        key_str = re.sub('g__','',key_raw)
        if key_str not in target_genus_dic:
            target_genus_dic[key_str] = [MAX,MIN]
        else:
            print("重复")
        

# 读取目标列信息

genus_lst = [] # 备份表头列
sample_name = ""
all_value_dic = {} # 用于存储
df_index = ['MAX', 'MIN']
with open('level-6-relative.csv',mode='rt',encoding='utf-8') as fh:
    '''
    index,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides,
    '''
    index_list  = [] # 目标的列索引,对应record
    for line in fh:
        record = re.split(',',line.strip())
        if line.startswith('index'):
            for index in range(1,len(record)-2):
                tmp = re.split(';',record[index])[-1]
                genus_name = re.sub('g__','',tmp)
                genus_lst.append(genus_name)
                if genus_name in target_genus_dic:
                    #print(genus_name)
                    #print(index)
                    index_list.append(index)
        else:
            sample_name = record[0]
            df_index.append(sample_name)
            for index in index_list:
                #print(index)
                genus_name = genus_lst[index]
                value = record[index]
                all_value_dic[genus_name]=value


for genus_name in target_genus_dic:
    if  genus_name in all_value_dic:
        value = all_value_dic[genus_name]
        target_genus_dic[genus_name].append(value)
    else:
        target_genus_dic[genus_name].append(MIN)

# print(target_genus_dic)
# print(list(genus_dic.values())[0])
# print(len(list(genus_dic.values())[0]))
df = pd.DataFrame(target_genus_dic)
df.index = df_index
df = df
print(df)
df.to_csv('xx.tsv',sep='\t')