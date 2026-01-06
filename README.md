# DNA Analysis with Claude

Analyze your genetic data through conversation with Claude. Works with 23andMe, AncestryDNA, MyHeritage, Nebula, and other DNA testing services. Run scripts, explore your genome, and get personalized explanations — all in natural dialogue.

## Why Claude + DNA?

Traditional DNA analysis tools give you static reports. With Claude, you get:

- **Interactive exploration** — Ask follow-up questions about any result
- **Plain English explanations** — No need to Google every SNP
- **Cross-category insights** — Connect findings across health, nutrition, ancestry
- **Personalized context** — Discuss how results apply to your specific situation
- **Research deep-dives** — Ask Claude to explain the science behind any marker

## Quick Start

1. Export raw data from your DNA testing service (see [Supported Services](#supported-dna-testing-services))
2. Place the file in `data/` folder
3. Open Claude Code in this project directory
4. Tell Claude: *"Update scripts to use my genome file: data/your_filename.txt"*
5. Start chatting:

```
You: Run the health analysis script

Claude: *runs health_analysis.py, shows results*

You: What does my MTHFR status mean for daily life?

Claude: *explains your specific genotype and practical implications*

You: Should I be concerned about the APOE result?

Claude: *provides context, explains penetrance, suggests next steps*
```

## Available Analysis Scripts

| Script | What it analyzes |
|--------|------------------|
| `health_analysis.py` | Cardiovascular, oncology, neurology, diabetes, pharmacogenomics |
| `ancestry_analysis.py` | Ethnic markers, haplogroups, pigmentation genes |
| `nutrition_analysis.py` | Lactose, gluten, caffeine, alcohol, vitamins metabolism |
| `carrier_status_analysis.py` | Hereditary condition carrier status |
| `cognitive_analysis.py` | Cognitive traits and predispositions |
| `longevity_analysis.py` | Longevity-associated variants |
| `psychology_analysis.py` | Psychological traits |
| `sports_fitness_analysis.py` | Athletic performance genetics |
| `sleep_chronotype_analysis.py` | Sleep patterns, circadian rhythm |
| `immunity_analysis.py` | Immune system markers |
| `detoxification_analysis.py` | Detox pathways, enzyme activity |
| `skin_analysis.py` | Skin characteristics |
| `vision_hearing_analysis.py` | Sensory genetics |
| `pain_sensitivity_analysis.py` | Pain perception |
| `reproductive_analysis.py` | Reproductive health |
| `physical_traits_analysis.py` | Physical characteristics |

## Example Conversations

### Exploring Health Results

```
You: Run health analysis and tell me the most important findings

You: I see I have a variant in CYP2D6. What medications should I be careful with?

You: My doctor wants to prescribe codeine. Is that safe for my genotype?
```

### Understanding Ancestry

```
You: Analyze my ancestry markers

You: What does the EDAR variant tell us about my ancestry?

You: How do my haplogroups connect to migration patterns?
```

### Optimizing Nutrition

```
You: What does my genome say about coffee metabolism?

You: Based on my variants, what supplements might actually help me?

You: I'm lactose intolerant symptoms-wise. What does my genetics show?
```

### Deep Dives

```
You: Explain the APOE gene and why it matters

You: What's the current research on rs1801133?

You: How reliable are these SNP associations? What's the evidence quality?
```

## Project Structure

```
dna-claude-analysis/
├── data/
│   ├── README.md            # Instructions for data placement
│   └── your_genome.txt      # ← Put your DNA file here
├── scripts/
│   └── *.py                 # Analysis scripts
├── reports/                 # Generated reports (created automatically)
├── .gitignore               # Keeps your data private
└── README.md
```

## Setup: Configuring Your Data File

**Important:** Before running scripts, you need to set the path to your genome file.

### Option 1: Ask Claude (Recommended)

Just tell Claude:
```
You: Update all scripts to use my genome file: data/my_dna_file.txt
```

Claude will update the `GENOME_FILE` variable in all scripts automatically.

### Option 2: Manual Setup

Each script has a `GENOME_FILE` variable at the top:

```python
# In scripts/health_analysis.py (and other scripts)
GENOME_FILE = f"{BASE_PATH}/data/your_genome_file.txt"  # ← Change this
```

Replace with your actual filename:
```python
GENOME_FILE = f"{BASE_PATH}/data/genome_John_Doe_v5_Full_20240101.txt"
```

## How It Works

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  Your Genome    │ ──▶ │  Python Script  │ ──▶ │  Markdown Report│
│  (600K SNPs)    │     │  (analyze SNPs) │     │  (findings)     │
└─────────────────┘     └─────────────────┘     └─────────────────┘
                                                        │
                                                        ▼
                                               ┌─────────────────┐
                                               │     Claude      │
                                               │  (explain,      │
                                               │   contextualize,│
                                               │   answer Qs)    │
                                               └─────────────────┘
```

## Tips for Best Results

1. **Start broad, then narrow** — Run a full analysis first, then ask about specific findings
2. **Ask "so what?"** — Claude can explain practical implications, not just raw data
3. **Request comparisons** — "How does my result compare to population averages?"
4. **Explore connections** — "Could my sleep issues be related to any of these variants?"
5. **Verify important findings** — Ask Claude about evidence quality and study sizes

## Important Notes

- **Not medical advice** — This is educational exploration, not diagnosis
- **Genetics ≠ destiny** — Environment and lifestyle matter more for most traits
- **Consult professionals** — For medical decisions, see a genetic counselor
- **Privacy first** — Your genome data stays local on your machine

## Requirements

- Python 3.6+
- Claude Code (or any Claude interface with file access)
- Raw DNA data from any major testing service

## Supported DNA Testing Services

Scripts work with raw data from multiple providers. Most use similar formats (rsID, chromosome, position, genotype).

### Consumer Genotyping (~600K-700K SNPs)

| Service | Coverage | Raw Data Export | Notes |
|---------|----------|-----------------|-------|
| [23andMe](https://www.23andme.com) | ~640K SNPs | Settings → 23andMe Data → Download | Most popular, good health reports |
| [AncestryDNA](https://www.ancestry.com/dna) | ~700K SNPs | Settings → Download DNA Data | Best for genealogy |
| [MyHeritage DNA](https://www.myheritage.com/dna) | ~700K SNPs | DNA → Manage DNA kits → Download | Good European coverage |
| [FamilyTreeDNA](https://www.familytreedna.com) | ~700K SNPs | myFTDNA → Download Raw Data | Y-DNA & mtDNA specialists |
| [Living DNA](https://www.livingdna.com) | ~650K SNPs | Account → Download | Best British Isles ancestry |

### Whole Genome Sequencing (WGS) — Full Coverage

| Service | Coverage | Price Range | Notes |
|---------|----------|-------------|-------|
| [Nebula Genomics](https://nebula.org) | 30x WGS | $250-300 | Privacy-focused, good reports |
| [Dante Labs](https://www.dantelabs.com) | 30x WGS | $200-400 | European company |
| [Sequencing.com](https://sequencing.com) | 30x WGS | $200-400 | Good app ecosystem |
| [Sano Genetics](https://sanogenetics.com) | 30x WGS | ~$300 | UK-based |

### Free / Research Options

| Service | Type | Notes |
|---------|------|-------|
| [All of Us](https://allofus.nih.gov) | WGS | NIH research program (US only) |
| [Our Future Health](https://ourfuturehealth.org.uk) | Genotyping | UK research program |
| [Genes for Good](https://genesforgood.sph.umich.edu) | Genotyping | University of Michigan (limited) |

### Where to Upload Existing Data

Already tested? Upload your raw data for additional analysis:

- [Promethease](https://promethease.com) — $12, comprehensive health report
- [GEDmatch](https://gedmatch.com) — Free, genealogy & ethnicity
- [Genetic Genie](https://geneticgenie.org) — Free, methylation & detox
- [Codegen.eu](https://codegen.eu) — Free, various health reports
- [DNA.Land](https://dna.land) — Free, ancestry & traits

## Getting Your Data

### 23andMe
1. Log into [23andMe](https://you.23andme.com)
2. Settings → 23andMe Data → Download Raw Data
3. Unzip and place in `data/` folder

### AncestryDNA
1. Log into [Ancestry](https://www.ancestry.com)
2. DNA → Settings → Download DNA Data
3. Confirm via email, download, place in `data/`

### MyHeritage
1. Log into [MyHeritage](https://www.myheritage.com)
2. DNA → Manage DNA kits → Actions → Download
3. Place in `data/` folder

### Other Services
Most services: Account Settings → Download Raw Data → Place in `data/`

## File Format

Scripts expect standard format (auto-detected):

```
# rsid          chromosome  position    genotype
rs4477212       1           82154       AA
rs3094315       1           752566      AG
```

Most services use this format. VCF files from WGS may need conversion (ask Claude for help).

---

*Explore your genome through conversation. Ask Claude anything.*
