

    ## Configure sra toolkit ##
# (implement a check for SRA toolkit)
# Download and extract NCBI SRA-toolkit: Ubuntu Lixux 64 bit archetecture version 2.11
#wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz;
#tar -xzvf sratoolkit.2.11.0-ubuntu64.tar.gz;
#rm sratoolkit.2.11.0-ubuntu64.tar.gz;

# Export to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

 ## configure (interactive) 
# opens, sets default settings and quits
echo fyoxyo | vdb-config -i;

# Update the prefetch download directory by editing SRA configuration file

echo '/repository/user/main/public/rt = '"\"$(pwd)/sra\"" >> $HOME/.ncbi/user-settings.mkfg;
echo '/repository/user/main/public/root = '"\"$(pwd)/sra\"" >> $HOME/.ncbi/user-settings.mkfg;


#cd sra
#vdb-config --prefetch-to-cwd;
#cd ..
