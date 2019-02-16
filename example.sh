# 1. Download an example input genome FASTA file
wget http://share.zhanxw.com/VAMPr/example/Klebsiella_pneumoniae_MGH_43.fasta

# 2. Download VAMPr protein database
perl VAMP_database.pl

# 3. Identify AMR genotypes
perl VAMP.pl Klebsiella_pneumoniae_MGH_43.fasta > Klebsiella_pneumoniae_MGH_43.VAMP.txt

# 4. Download VAMPr prediction model - you can find the download links in http://share.zhanxw.com/VAMPr/VAMPr.cgi?show=model
wget http://share.zhanxw.com/VAMPr/VAMP_model/573.amikacin.VAMP_model.RData
wget http://share.zhanxw.com/VAMPr/VAMP_model/573.cefepime.VAMP_model.RData

# 5. Predict AMR phenotypes
perl VAMP_predict.pl 573.amikacin.VAMP_model.RData amikacin=Klebsiella_pneumoniae_MGH_43.VAMP.txt
perl VAMP_predict.pl 573.cefepime.VAMP_model.RData cefepime=Klebsiella_pneumoniae_MGH_43.VAMP.txt
