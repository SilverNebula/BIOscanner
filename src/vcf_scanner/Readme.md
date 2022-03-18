## Vcf scanner tool
基于SnpSift工具将突变位点数据注释到ClinVar等数据库，并筛选结果

### requirements：
#### Java
```
(jdk >=1.8)
```
#### Python packages:
```
numpy
pyvcf3
pandas>='1.0.3'
pyhgvs
```

流程：
1.准备待注释的源vcf文件和所用数据库文件。
2.使用java工具SnpSift注释变异数据文件：

关于SnpSift的使用参考了：https://zhuanlan.zhihu.com/p/185970208

```
java -Xmx1g -jar ./SnpSift.jar annotate -v [clinvar database file] [source vcf file] > [result vcf file]
```

3.使用本工具vcf scanner扫描被注释后的vcf文件
```
python vcf_scanner.py [file name]
```