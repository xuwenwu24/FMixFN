#This is a reference text for file processing, we assume that the test file is a genotype binary file, 
#including the test.bed, test.bim, and test.fam files.First we get the clean file through QC. In this 
#example, the clean file is divided into training and validation sets, Note: These two files cannot contain 
#pure-fit bits, such as a column of data that cannot be 0, 1 or 2 after the file is converted to 0, 1, 2.
#The detailed processing format sits on each file is detailed in user manual.
#Data quality control
/lustre/softs/plink1.9 --bfile test --maf 0.05 --hwe 0.000001 --geno 0.1 --mind 0.1 --make-bed --out clean
##cross validation
cat clean.fam | sed -n '1,100p' | awk '{print $1,$2}' > validation
cat clean.fam | sed -n '101,$p' | awk '{print $1,$2}' > train
##converted data to 0,1,2 format
##added the column name to the data
/lustre/softs/plink1.9 --bfile clean --keep validation --recodeA --out validation_tmp
/lustre/softs/plink1.9 --bfile clean --keep train --recodeA --out train_tmp
##Turn training genotype files into 0, 1, 2 forms and add column names
sed -i '1d' train_tmp.raw  
awk '{printf("%d %s\n",NR,$0)}' train_tmp.raw > tmp
awk '{printf $1; for(i=8;i<=NF;i+=1) {printf " " $i}; printf "\n"}' tmp > training_genotype.txt
sed -n '1p' training_genotype.txt > head
cat training_genotype.txt >> head #add a column name
mv head training_genotype.txt 
##Turn validation genotype files into 0, 1, 2 forms and add column names
sed -i '1d' validation_tmp.raw
awk '{printf("%d %s\n",NR,$0)}' validation_tmp.raw > tmp
awk '{printf $1; for(i=8;i<=NF;i+=1) {printf " " $i}; printf "\n"}' tmp > validation_genotype.txt
sed -n '1p' validation_genotype.txt > head
cat validation_genotype.txt >> head 
mv head validation_genotype.txt 
##genotype files were ok, preparing the phenotype file
awk '{print $6,$5}' clean.fam > tmp
cat tmp | sed -n '1,100p' > validation_fix
awk '{printf("%d %s\n",NR,$0)}' validation_fix > validation_phenotype.txt
cat tmp | sed -n '101,$p'  > train_fix
awk '{printf("%d %s\n",NR,$0)}' train_fix > training_phenotype.txt
##phenotype is ok, preparing the parameter file 
cat train_tmp.raw | wc | awk '{print $1}' > test.txt #the number of training individuals
cat train_tmp.raw | awk '{print NF-6}' | sed -n '1p' >> test.txt #add the number of markers
echo 1 >> test.txt #add the number of fixed effects, Suppose the number of fixed effects is 1
echo 0.25 >> test.txt #add PI value
echo 0.61 >> test.txt #add the heritability value, Suppose the heritability value is 0.61
cat validation_tmp.raw | wc | awk '{print $1}'  >> test.txt #add the number of markers
echo 1 >> test.txt #Whether calculate the predictive ability (input 1 means to calculate, input 0 means not calculated)

##Start running FCF-MixP with all files ready
##Run
./FCF-MixP





