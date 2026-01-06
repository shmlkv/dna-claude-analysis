# Place Your DNA Data Here

Put your raw DNA file in this folder.

## Expected File Format

Text file with tab-separated values:

```
# rsid          chromosome  position    genotype
rs4477212       1           82154       AA
rs3094315       1           752566      AG
rs1801133       1           11856378    CT
```

## Supported File Types

- `.txt` — 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA raw exports
- `.csv` — Some services export as CSV
- `.vcf` — VCF files from whole genome sequencing (may need conversion)

## How to Get Your Data

### Consumer Genotyping Services

| Service | How to Download |
|---------|-----------------|
| **AncestryDNA** | ancestry.com → DNA → Settings → Download DNA Data |
| **23andMe** | you.23andme.com → Settings → 23andMe Data → Download |
| **MyHeritage** | myheritage.com → DNA → Manage DNA kits → Download |
| **FamilyTreeDNA** | familytreedna.com → myFTDNA → Download Raw Data |
| **Living DNA** | livingdna.com → Account → Download |

### Whole Genome Sequencing

| Service | How to Download |
|---------|-----------------|
| **Nebula Genomics** | nebula.org → Account → Download Raw Data |
| **Dante Labs** | dantelabs.com → My Account → Download |
| **Sequencing.com** | sequencing.com → My Files → Download |
| **Sano Genetics** | sanogenetics.com → Dashboard → Download |

### Research Programs

| Program | Notes |
|---------|-------|
| **All of Us** (NIH) | allofus.nih.gov — US residents |
| **Our Future Health** | ourfuturehealth.org.uk — UK residents |
| **Genes for Good** | University of Michigan program |

### Converting Other Formats

Have a VCF file from WGS? Ask Claude:
```
You: Convert my VCF file to 23andMe format
```

## After Placing Your File

Tell Claude to configure the scripts:
```
You: Update all scripts to use my genome file: data/my_filename.txt
```

## Privacy & Security

- This folder is in `.gitignore` — your data won't be committed to git
- Keep your DNA data private and secure
- Never share raw DNA files publicly
- Consider encrypting backups of your genetic data
