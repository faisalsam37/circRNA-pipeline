# circRNA Pipeline Setup

This pipeline offers a guided, interactive experience to help users easily process RNA-seq data for circular RNA (circRNA) detection.

Once the `hc_circ_rna_hardcoded.sh` script is launched, the user is prompted to:

- Choose one of the four supported circRNA detection tools:
  - find_circ
  - CIRI2
  - CLEAR / CIRCexplorer3
  - circRNA_finder
- Enter an SRA accession number corresponding to the RNA-seq dataset to be processed.

The script will then automatically download the dataset, perform preprocessing, align reads, and run the selected circRNA detection tool.

## STEP1: Make sure target directory exists

```bash
cd /mnt/Data/research
# Ensure environment.yml is present
ls environment.yml
```


## STEP2: Make sure target directory exists

```bash
cd /mnt/Data/research

# Ensure environment.yml is present
ls environment.yml

```bash
cd /mnt/Data/research
## Ensure environment.yml is present
ls environment.yml
```

## STEP3: Create and activate conda environment
```bash
conda env create -f environment.yml
conda activate circRNA_pipeline
```

# Create and navigate to the tools directory
```bash
mkdir -p /mnt/Data/research/tools
cd /mnt/Data/research/tools

### CLEAR / CIRCexplorer3
git clone https://github.com/YangLab/CLEAR.git
cd CLEAR
python ./setup.py install
cd ..

### find_circ
git clone https://github.com/marvin-jens/find_circ.git
cd find_circ
# NOTE: No install needed; it's a script-based tool
cd ..

### CIRI2
wget https://downloads.sourceforge.net/project/ciri/CIRI2/CIRI_v2.0.6.zip
unzip CIRI_v2.0.6.zip -d CIRI2
cd CIRI2
# CIRI2 is a Perl script — no install needed
cd ..

### circRNA_finder
conda install -c bioconda circrna_finder
which circrna_finder
cp $(which circrna_finder) /mnt/Data/research/tools/
```
