
import re
import os
from pathlib  import Path
import pandas as pd
import sys

'''
批量读取，删除固定的列

筛选Significant 的内容 提取基因列
合并为一个dataframe

df = pd.DataFrame({'col1': list1, 'col2': list2})
print(df)
写文件

'''

def GetAllFilePaths(pwd,wildcard='*'):
    '''
    获取目录下文件全路径，通配符检索特定文件名，返回列表
    param: str  "pwd"
    return:dirname pathlab_obj
    return:list [ str ]
    #https://zhuanlan.zhihu.com/p/36711862
    #https://www.cnblogs.com/sigai/p/8074329.html
    '''
    files_lst = []
    target_path=Path(pwd)
    for child in target_path.rglob(wildcard):
        if child.is_dir():
            pass
        elif child.is_file():
            files_lst.append(child)
    return files_lst

def nan_fill(dic):
    '''
    dic= {
        'a':[1,2,3],
        'b':[5,],
        'c':[4,6,7,8]
        }
    '''
    max_len = max([len(v) for v in dic.values()])
    for k in dic:
        lst_len = len(dic[k])
        if len(dic[k]) < max_len:
            dic[k] += [pd.NA,]*(max_len-lst_len)
    return dic


def get_union_gene(sig_gene_dic):
    union_gene_dic = {}
    other_set_gene_dic = {}
    sample_list = list(sig_gene_dic.keys())
    for sample in sample_list:
        if sample not in union_gene_dic:
            # init variable
            other_set_gene_dic[sample] = set()
            union_gene_dic[sample] = set()
        for other_sample in sig_gene_dic:
            if other_sample != sample:
                other_set_gene_dic[sample] = other_set_gene_dic[sample].union(sig_gene_dic[other_sample])

    for  sample in sample_list:
        set1 = set(sig_gene_dic[sample])
        other_set2 = other_set_gene_dic[sample]
        union_gene_dic[sample] = list(set1.difference(other_set2))

    return union_gene_dic

def main():

    pwd = Path.cwd()
    os.chdir(pwd)

    
    file_list = GetAllFilePaths(pwd,wildcard='*.anno.xls')
    sig_gene_dic = {} # 存储每一个样品中显著差异的基因名list
    dataframe_dic = {} # sample : dataframe
    sample_name_lst = []
    all_gene_set = set() # 用于输出长度一致的表格
    for file_path in file_list:
        sample_name = re.sub(".anno","", file_path.stem)
        sample_name_lst.append(sample_name)
        df = pd.read_csv(file_path,sep='\t',index_col=None)
        #print(df)
        # 获取名字
        sig_df = df[df['Significant'] == 'yes']
        #print(sig_df)
        if sample_name not in sig_gene_dic:
            # 获得所有样品的基因并集，获得并集
            all_gene_set = all_gene_set.union(set(sig_df['Gene'].to_list()))
            # 获得每个样品的显著基因
            sig_gene_dic[sample_name] = set(sig_df['Gene'].to_list())
        else:
            print("重复样品：",sample_name)
        # 存储去掉前几列的数据，存储
        new_df = sig_df.drop(df.columns[1:14], axis=1, inplace=False)
        if sample_name not in dataframe_dic:
            dataframe_dic[sample_name] = new_df
        else:
            print("重复样品：",sample_name)
    
    # print(all_gene_set)
    # print(len(all_gene_set))
    # 计算样品的补集，独有的
    '''
    需要使用 另外3个样品的基因集合与目前的集合计算 补集
    
    '''
    union_gene_dic = get_union_gene(sig_gene_dic) # sample_name: union_gene_lst
    '''
    全部基因则是 所有显著性基因的并集
    '''
    all_gene_lst = list(all_gene_set)
    #print(all_gene_set)
    venn_bool_dic = {'value':[x for x in range(1,len(all_gene_lst)+1)],}
    # 对每一个datafame 筛选基因列，并分别命名输出
    for sample_name in sample_name_lst:
        if sample_name not in venn_bool_dic:
            venn_bool_dic[sample_name] = ['FALSE',]*len(all_gene_lst)
        #print(venn_bool_dic[sample_name])
        df = dataframe_dic[sample_name]
        gene_tmp_lst = df['Gene'].to_list()
        for index in range(len(all_gene_lst)):
            #print(all_gene_lst[index])
            if all_gene_lst[index] in gene_tmp_lst:
                venn_bool_dic[sample_name][index] = "TRUE"
        #print(union_gene_dic[sample_name])
        out_file_name = pwd.joinpath(sample_name + '.anno.Union.xls')
        new_df = df[df['Gene'].isin(union_gene_dic[sample_name])]
        new_df.to_csv(out_file_name,sep='\t',encoding="utf-8", index=False, header=True)

    '''
    Element  name1 name2 name3 name4
    '''
    #print(venn_bool_dic)
    # 合并union_gene_dic 转为dataframe 输出
    new_data_dic = nan_fill(venn_bool_dic)
    venn_plot_df = pd.DataFrame(new_data_dic)
    #print(venn_bool_dic)
    #venn_plot_df = pd.DataFrame(pd.DataFrame.from_dict(union_gene_dic, orient='index').values.T, columns=list(union_gene_dic.keys()))
    venn_plot_df.to_csv(pwd.joinpath('venn_plot.tsv'),sep='\t',encoding="utf-8", index=False, header=True)



if __name__ == '__main__':
    main()