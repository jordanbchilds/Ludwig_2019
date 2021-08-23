#!/bin/bash


    ## Download and configure sra toolkit ##

  # MANUAL CONFIGURATION NECESSARY
  # EXECUTE THE FOLLOWING COMMANDS IN YOUR TERMINAL
  # When `vdb-config -i` opens an interactive configuration window, follow the instructions to exit, and continue to execute the following lines:


# Download and extract NCBI SRA-toolkit: Ubuntu Lixux 64 bit archetecture version 2.11
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz;
tar -xzvf sratoolkit.2.11.0-ubuntu64.tar.gz;
rm sratoolkit.2.11.0-ubuntu64.tar.gz;

# Export to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

 ## configure (interactive) 
# type "fyoxyo" when interactive display opens
 vdb-config -i;

# Update the prefetch download directory by editing SRA configuration file

echo '/repository/user/main/public/rt = '"\"$(pwd)/sra\"" >> $HOME/.ncbi/user-settings.mkfg;
echo '/repository/user/main/public/root = '"\"$(pwd)/sra\"" >> $HOME/.ncbi/user-settings.mkfg;


