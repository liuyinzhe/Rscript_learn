import os,re
from pathlib import Path
import pandas as pd

'''
输入：
    rno00001.keg
    MK.differential_gene.tsv
    M_vs_K.KEGG_enrich.tsv

输出：
    MK.KEGG_bar_plot.tsv
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


def main():
    name='ko'
    ## 跳转脚本所在目录
    #pwd = os.path.split(os.path.realpath(__file__))[0]
    #print(pwd)
    #pwd = os.getcwd()
    #pwd = Path(pwd)
    pwd = Path.cwd()
    os.chdir(pwd)
    
    LV2_dic = {}
    ko_dir = Path("/data/home/liuyinzhe/project/RNA_seq/database")
    ko_file = ko_dir.joinpath("hsa00001.keg")
    with open(ko_file,mode='rt',encoding='utf-8') as fh:
        class_type = ""
        for line in fh:
            #print(line)
            if line.startswith("D") or line.startswith("#") or line.startswith("%") or line.startswith("!"):
                continue
            #print('xxx')
            if line.startswith("A") and len(line)>2:
                record = re.split("\s+?",line.strip(),maxsplit=1)
                class_type = record[1]
                #print(class_type)
                continue
            if "[PATH:" in line:
                #print(line)
                match_obj = re.search("\[PATH\:(\w{3,})\]",line.strip())
                result_str = match_obj.group(1)
                if result_str not in LV2_dic:
                    LV2_dic[result_str] = class_type
                # if 'rno04020' in line:
                #     print(line)
                #     print("#"+result_str+"#")
                #     print(result_str,"#",class_type)
                #     print(LV2_dic[result_str])
    #print(LV2_dic)
            
    top_num = 9999 #*2 
    pvalue_cutoff = 1 #0.05
    all_KEGG_info = []
    '''
    description(LV2),significance,class(LV1),gene_count,pvalue
    '''
    Down_kegg_dir = pwd.joinpath("KEGG","Down")
    Down_kegg_file_lst = GetAllFilePaths(Down_kegg_dir,wildcard="*.Down_KEGG_enrich.tsv")
    
    down_plot_file = pwd.joinpath("KEGG","Down","KEGG_bar_plot.tsv")

    Up_kegg_dir= pwd.joinpath("KEGG","Up")
    Up_kegg_file_lst = GetAllFilePaths(Up_kegg_dir,wildcard="*.Up_KEGG_enrich.tsv")
    
    up_plot_file = pwd.joinpath("KEGG","Up","KEGG_bar_plot.tsv")
    if len(Down_kegg_file_lst)>0:
        Down_kegg_file = Down_kegg_file_lst[0]
        with open(Down_kegg_file,mode='rt',encoding='utf-8') as fh:
            for raw_line in fh:
                line= re.sub('\"','',raw_line.strip())
                if line.startswith("ID"):
                        continue
                record = re.split('\t',line)
                '''
                "ID"	"Description"	"GeneRatio"	"BgRatio"	"pvalue"	"p.adjust"	"qvalue"	"geneID"	"Count"
                0            1           2              3          4          5            6          7           8 
                '''
                description = record[1] # 当前级别注释描述
                gene_ratio = record[2]
                bg_ratio =  record[3]
                pvalue = float(record[4]) # pvalue
                qvalue = float(record[5]) # p.adjust
                #geneid = record[7] # gene1/gene2
                #geneid_lst = re.split("/",geneid)
                kegg_ID = record[0]
                #print(kegg_ID)
                if kegg_ID not in LV2_dic:
                    #print(kegg_ID)
                    continue
                kegg_class = LV2_dic[kegg_ID] # level 2 分类描述

                gene_count = int(record[8])
                # pvalue 阈值筛选
                if pvalue>pvalue_cutoff:
                    continue
                all_KEGG_info.append([description,kegg_class,gene_ratio,bg_ratio,pvalue,qvalue,gene_count])
                
            # pvalue,qvalue 降序排序,gene_count 升序排序
            all_KEGG_sorted = sorted(all_KEGG_info,key=lambda x:(x[4],x[5],-x[6]))

        count = 0
        with open(down_plot_file,mode='wt',encoding='utf-8') as out:
            out.write('description\tkegg_class\tGeneRatio\tBgRatio\tp_value\tq_value\tgene_count\n')
            for x in all_KEGG_sorted:
                count +=1
                #print(count,x)
                out.write('\t'.join(list(map(str,x)))+'\n')
                if count == top_num:
                    break
    # UP
    if len(Up_kegg_file_lst)>0:
        Up_kegg_file = Up_kegg_file_lst[0]
        with open(Up_kegg_file,mode='rt',encoding='utf-8') as fh:
            for raw_line in fh:
                line= re.sub('\"','',raw_line.strip())
                if line.startswith("ID"):
                        continue
                record = re.split('\t',line)
                '''
                "ID"	"Description"	"GeneRatio"	"BgRatio"	"pvalue"	"p.adjust"	"qvalue"	"geneID"	"Count"
                0            1           2              3          4          5            6          7           8 
                '''
                description = record[1] # 当前级别注释描述
                gene_ratio = record[2]
                bg_ratio =  record[3]
                pvalue = float(record[4]) # pvalue
                qvalue = float(record[5]) # p.adjust
                #geneid = record[7] # gene1/gene2
                #geneid_lst = re.split("/",geneid)
                kegg_ID = record[0]
                #print(kegg_ID)
                if kegg_ID not in LV2_dic:
                    #print(kegg_ID)
                    continue
                kegg_class = LV2_dic[kegg_ID] # level 2 分类描述

                gene_count = int(record[8])
                # pvalue 阈值筛选
                if pvalue>pvalue_cutoff:
                    continue
                all_KEGG_info.append([description,kegg_class,gene_ratio,bg_ratio,pvalue,qvalue,gene_count])
                
            # pvalue,qvalue 降序排序,gene_count 升序排序
            all_KEGG_sorted = sorted(all_KEGG_info,key=lambda x:(x[4],x[5],-x[6]))

        count = 0
        with open(up_plot_file,mode='wt',encoding='utf-8') as out:
            out.write('description\tkegg_class\tGeneRatio\tBgRatio\tp_value\tq_value\tgene_count\n')
            for x in all_KEGG_sorted:
                count +=1
                #print(count,x)
                out.write('\t'.join(list(map(str,x)))+'\n')
                if count == top_num:
                    break
    
        
if __name__ == '__main__':
    main()
