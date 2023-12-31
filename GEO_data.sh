#Create data folder
mkdir data

#Download file from NCBI server into data folder
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72857/suppl/GSE72857%5Fumitab%2Etxt%2Egz -o data/GSE72857_umitab.txt.gz

#Decompress it
gunzip data/GSE72857_umitab.txt.gz
