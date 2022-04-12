
cd /media/datn/data/LmTag_revise

git clone https://github.com/datngu/LmTag_dev.git

cd LmTag_dev
# build the C++ program
make
# export to PATH
export PATH=$PWD/bin:$PATH
cd ..

##
# download dataset

mkdir EAS_running_time

cd EAS_running_time

wget https://zenodo.org/record/5807198/files/chr10_EAS.vcf.gz?download=1 -O chr10_EAS.vcf.gz

wget https://zenodo.org/record/5807198/files/chr10_EAS_CADD.txt?download=1 -O chr10_EAS_CADD.txt

wget https://zenodo.org/record/5807198/files/VIP_GWAS_CLINVAR_ALL.txt?download=1 -O VIP_GWAS_CLINVAR_ALL.txt

wget https://zenodo.org/record/5807198/files/chr10_EAS.tar.gz?download=1 -O lmtag_test/chr10_EAS.tar.gz

tar -xvzf chr10_EAS.tar.gz


# build model

# assumming you are running unix and $PWD/lmtag_test: is the path to lmtag_test directory that is mounted to docker containter as /lmtag_test:

model_pipeline.sh \
      -v chr10_EAS_processed.vcf.gz \
      -r chr10_EAS_hg38_high_cov \
      -n 32970 \
      -p 16 \
      -o LmTag_model_test.Rdata


# running time K values


for K in 1 5 10 20 30 50 100 200 500 1000 1500 2000
do
	LmTag_pipeline.sh -v chr10_EAS_processed.vcf.gz -m LmTag_model_test.Rdata -s chr10_EAS_CADD.txt -V VIP_GWAS_CLINVAR_ALL.txt -k $K -M linear -o test_LmTag_model_linear_${K}.txt
done






###############################################
################################################
#################################################

### Installing and running other tools


cd /media/datn/data/LmTag_revise

git clone https://github.com/datngu/TagSNP_evaluation

# Installing tools

# TagIt
cd TagSNP_evaluation/TagIt
bash configure.sh
TagIt_tool=$PWD/TagIt_pipeline.sh
cd ../..

# FastTagger

cd TagSNP_evaluation/FastTagger 
bash configure.sh
FastTagger_tool=$PWD/FastTagger_pipeline.sh
cd ../..

# EQ_uniform

EQ_uniform_tool=$PWD/TagSNP_evaluation/EQ_uniform/EQ_uniform.R

# EQ_MAF

EQ_MAF_tool=$PWD/TagSNP_evaluation/EQ_MAF/EQ_MAF.R

### Running tools
cd EAS_running_time

# TagIt
start_time=$SECONDS

bash $TagIt_tool chr10_EAS_processed.vcf.gz chr10_EAS_TagIt
head -n 32000 chr10_EAS_TagIt/chr10_EAS_TagIt_tags_cleaned.txt > chr10_EAS_TagIt_32000_cleaned.txt

elapsed=$(( SECONDS - start_time ))
eval "echo Running time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')" > TagIt.time.log
echo DONE!


# FastTagger
cd EAS_running_time

start_time=$SECONDS

bash $FastTagger_tool chr10_EAS_processed.vcf.gz 32000 chr10_EAS_FastTagger

elapsed=$(( SECONDS - start_time ))
eval "echo Running time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')" > FastTagger.time.log
echo DONE!

# EQ_uniform
cd EAS_running_time

start_time=$SECONDS

Rscript $EQ_uniform_tool vcf=chr10_EAS_processed.vcf.gz size= 32000 out=chr10_EAS_32000_EQ_uniform.txt
elapsed=$(( SECONDS - start_time ))
eval "echo Running time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')" > EQ_uniform.time.log
echo DONE!

# EQ_uniform
cd EAS_running_time

start_time=$SECONDS

Rscript $EQ_MAF_tool vcf=chr10_EAS_processed.vcf.gz size=32000 out=chr10_EAS_32000_EQ_MAF.txt
elapsed=$(( SECONDS - start_time ))
eval "echo Running time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')" > EQ_MAF.time.log
echo DONE!






















