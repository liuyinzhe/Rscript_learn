import pandas as pd

'''
*_multianno.xls 文件转为 annovarToMaf 输入的 annovar 文件格式
# 需要注意 注释时用的是 ensGene 还是 refGene ,annovarToMaf() 函数参数需要调整 table = 'ensGene' 或者 table = 'refGene'
'''

usecols_lst = [
    'Chr',
    'Start',
    'End',
    'Ref',
    'Alt',
    'Gene.ensGene', # 'Gene.refGene',
    'GeneDetail.ensGene', # 'GeneDetail.refGene',
    'ExonicFunc.ensGene', # 'ExonicFunc.refGene',
    'AAChange.ensGene', # 'AAChange.refGene',
    #'Tumor_Sample_Barcode',
    'Func.ensGene' # 'Func.refGene'
]
sample_name ='A'

multianno_snp_df = pd.read_csv(sample_name+'.SNP.hg19_multianno.xls',sep='\t',encoding='utf-8',usecols=usecols_lst)
# print(multianno_snp_df)
multianno_snp_df['Tumor_Sample_Barcode'] = sample_name


multianno_indel_df = pd.read_csv(sample_name+'.INDEL.hg19_multianno.xls',sep='\t',encoding='utf-8',usecols=usecols_lst)
# print(multianno_indel_df)
multianno_indel_df['Tumor_Sample_Barcode'] = sample_name

multianno_df = pd.concat([multianno_snp_df, multianno_indel_df], ignore_index = True, join='outer') #处理其他轴上的索引：默认为'outer'，取并集；当join='inner'时，取交集
multianno_df.to_csv(sample_name+'.annovar',sep='\t',index=False)

'''
Pandas中的连接函数汇总
https://blog.csdn.net/Hour__/article/details/118682204
'''
