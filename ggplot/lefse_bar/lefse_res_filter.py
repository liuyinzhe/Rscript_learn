import re

# Biomarker_names Logarithm_value Groups LDA_value(log10) P_value 
# 0                 1               2       3       4
line_lst = []
with open("lefse.res",mode='rt',encoding='utf-8') as fh:
    for line in fh:
        record = re.split("\t",line.strip())

        if len(record) <=3:
            continue
        group = record[2]
        if not group:
            continue
        LDA_score = record[3] # 用于限制 LDA_score 大小，软件绘制时默认只读取小于2的进行绘图
        #print(LDA_score)
        if float(LDA_score) < 2:
            continue
        record[0] = re.split('\.',record[0])[-1]
        # 物种层级筛选
        # # 界 Kingdom  门 Phylum 纲 Class 目 Order 科 Family Genus 种 Species
        # if not record[0].startswith(s__):
        #    continue
        line_lst.append(record)


with open("lefse.res.tsv",mode='wt',encoding='utf-8') as out:
    count=0
    out.write('\t'.join(['biomarker','logarithm','groups','LDA_score','p_value'])+'\n')
    for lst in sorted(line_lst,key=lambda x:x[3],reverse = True):
        out.write('\t'.join(lst)+'\n')
        # 总数限制
        # count+=1
        # if count ==20 :
        #     break

positive_tag = 'QM'
negative_tag = 'SN'

with open("lefse.res.plot.tsv",mode='wt',encoding='utf-8') as out:
    count=0
    out.write('\t'.join(['biomarker','logarithm','groups','LDA_score','p_value'])+'\n')
    for lst in sorted(line_lst,key=lambda x:x[3],reverse = True):
        group = lst[2]
        if group == negative_tag:
            lst[3] = str(float(lst[3]) * -1)
        out.write('\t'.join(lst)+'\n')
        # 总数限制
        # count+=1
        # if count ==20 :
        #     break
