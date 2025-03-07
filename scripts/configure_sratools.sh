#!/bin/bash

# SRAtoolkit must be configured interactively: when installing on rocket do not submit this as a batch job, but run on a login terminal or interactive node (srun)
# To do this, execute this file **configure_sratools.sh** line by line from the login node terminal (ie. do not submit script to SLURM) by reading, then pasting and executing all commands below. When prompted interactively set default configuration by inputting: "f","y","o","x","y","o".
 

    ## Download and configure sra toolkit ##

  # MANUAL CONFIGURATION NECESSARY
  # EXECUTE THE FOLLOWING COMMANDS LINE BY LINE IN THE TERMINAL FROM THE PROJECT ROOT DIR "Ludwig_2019/"
  
# MAKE SURE to execute from the root dir of the project "Ludwig_2019/", (not "Ludwig_2019/scripts/")


# Download and extract NCBI SRA-toolkit: Ubuntu Lixux 64 bit archetecture version 2.11
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz;
tar -xzvf sratoolkit.2.11.0-ubuntu64.tar.gz;
rm sratoolkit.2.11.0-ubuntu64.tar.gz;

# Export to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

  ## Configure (interactive) ## 

# When `vdb-config -i` opens an interactive configuration window, follow the instructions to set the defaults and exit, and continue line by line execution of this script.
# Execute line and follow the instructions to set the defaults and exit (input "fyoxyo")
vdb-config -i;
# This will create an SRA configuration file in $HOME/.ncbi/

# Update the prefetch download directory by editing SRA configuration file. (~/.ncbi/) 
echo '/repository/user/main/public/rt = '"\"$(pwd)/sra\"" >> $HOME/.ncbi/user-settings.mkfg;
echo '/repository/user/main/public/root = '"\"$(pwd)/sra\"" >> $HOME/.ncbi/user-settings.mkfg;

