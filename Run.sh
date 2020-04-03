#This is a example for FCF-MixP
#The files to be prepared include: training_genotype.txt,training_phenotype.txt,validation_genotype.txt,validation_phenotype.txt and parameter.txt.
#In this example, there are 739 individuals and 33901 maekers in the training population, and the validation population has 100 individuals. The heritability of the trait is 0.3.
##Run
gfortran -o FCF-MixP target.f90 -ffree-line-length-none -m32
./FCF-MixP

