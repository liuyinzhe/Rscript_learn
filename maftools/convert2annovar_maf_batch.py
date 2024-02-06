import pandas as pd
import re
'''
*_multianno.txt 文件转为 annovarToMaf 输入的 annovar 文件格式
# 需要注意 注释时用的是 ensGene 还是 refGene ,annovarToMaf() 函数参数需要调整 table = 'ensGene' 或者 table = 'refGene'
'''

def convert2annovar(sample_name,new_name,ref_ver='hg19'):

    usecols_lst = [
        'Chr',
        'Start',
        'End',
        'Ref',
        'Alt',
        'Gene.refGene',
        'GeneDetail.refGene',
        'ExonicFunc.refGene',
        'AAChange.refGene',
        'Func.refGene'
        # 'Gene.ensGene', # 'Gene.refGene',
        # 'GeneDetail.ensGene', # 'GeneDetail.refGene',
        # 'ExonicFunc.ensGene', # 'ExonicFunc.refGene',
        # 'AAChange.ensGene', # 'AAChange.refGene',
        #'Tumor_Sample_Barcode',
        #'Func.ensGene' # 'Func.refGene'
    ]
    try:
        multianno_snp_df = pd.read_csv(sample_name+'.snp.'+ref_ver+'_multianno.txt',sep='\t',encoding='utf-8',usecols=usecols_lst)
        #print(multianno_snp_df)
        multianno_snp_df['Tumor_Sample_Barcode'] = new_name
    except pd.errors.EmptyDataError:
        multianno_snp_df = pd.DataFrame()

    try:
        multianno_indel_df = pd.read_csv(sample_name+'.indel.'+ref_ver+'_multianno.txt',sep='\t',encoding='utf-8',usecols=usecols_lst)
        #print(multianno_indel_df)
        multianno_indel_df['Tumor_Sample_Barcode'] = new_name
    except pd.errors.EmptyDataError:
        multianno_indel_df = pd.DataFrame()

    return multianno_snp_df,multianno_indel_df

def main():
    all_annovar_df_list = []
    name_dic = {}
    with open('all.pair_info.txt',mode='rt',encoding='utf-8') as fh,open('clinical.txt',mode='wt',encoding='utf-8') as out :
        out.write('Tumor_Sample_Barcode\tSample_organization\n')
        for line in fh:
            '''
            uterus_G8_P1	PC_24_PC_33

            Tumor_Sample_Barcode	Sample_organization
            '''
            record = re.split('\t',line.strip())
            sample_name = record[0]
            sample_id = record[1]
            organization = re.split('\t',sample_name)[0]
            name_dic[sample_name] = sample_id
            out.write('\t'.join([sample_name,organization])+'\n')
            snp_df,indel_df = convert2annovar(sample_id,sample_name,ref_ver='hg19')

            if snp_df.shape != (0,0) and indel_df.shape != (0,0) :
                multianno_df = pd.concat([snp_df,indel_df], ignore_index = True, join='outer')
                multianno_df.to_csv(sample_name+'.annovar',sep='\t',index=False)
            elif snp_df.shape == (0,0) and indel_df.shape == (0,0) :
                continue
            elif indel_df.shape == (0,0) and snp_df.shape != (0,0) :
                multianno_df = snp_df
                multianno_df.to_csv(sample_name+'.annovar',sep='\t',index=False)
            elif indel_df.shape != (0,0) and snp_df.shape == (0,0) :
                multianno_df = indel_df
                multianno_df.to_csv(sample_name+'.annovar',sep='\t',index=False)
            all_annovar_df_list.append(multianno_df)
    
    multianno_all = pd.concat(all_annovar_df_list, ignore_index = True, join='outer') 
    multianno_all.to_csv('all.annovar',sep='\t',index=False)

if __name__ == '__main__':
    main()
