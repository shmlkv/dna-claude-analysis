# DNA Analysis Project — Claude Instructions

## Project Overview

Personal genome analysis toolkit. Scripts analyze raw DNA data and output markdown reports. The final step is generating a single-page HTML visualization in terminal/hacker aesthetic.

## Workflow

```
1. User places DNA file in data/
2. Configure scripts: GENOME_FILE variable
3. Run analysis scripts → markdown reports in reports/
4. Generate webpage from reports → webpage/dna_terminal.html
```

## Key Commands

### Run All Analysis
```bash
cd scripts && for f in *_analysis.py; do python "$f"; done
```

### Generate Webpage
After running scripts, user asks:
```
"Generate DNA Terminal webpage from the reports"
```

## Webpage Generation

Reference `webpage/STYLE_GUIDE.md` for exact styling. Key points:

- **Single HTML file** with inline CSS/JS
- **Terminal aesthetic**: JetBrains Mono, green-on-black, Matrix vibes
- **Fixed header** with navigation
- **Sections** for each analysis category
- **Color coding**: green=good, amber=warning, red=risk

### Section IDs
| Section | ID |
|---------|-----|
| Status Overview | status |
| Ancestry | ancestry |
| Health Risks | health |
| Nutrition | nutrition |
| Sports/Fitness | sports |
| Psychology | psychology |
| Cognitive | cognitive |
| Longevity | longevity |
| Sleep | sleep |
| Immunity | immunity |
| Pain Sensitivity | pain |
| Detoxification | detox |
| Skin | skin |
| Vision/Hearing | vision |
| Physical Traits | traits |
| Pharmacogenomics | pharma |
| Carrier Status | carrier |

### Data Flow
```
reports/*.md → Parse findings → HTML sections → webpage/dna_terminal.html
```

## File Locations

| Path | Purpose |
|------|---------|
| `data/` | Raw DNA files (gitignored) |
| `scripts/` | Python analysis scripts |
| `reports/` | Generated markdown reports |
| `webpage/` | Style guide + generated HTML |

## Important Rules

1. **Privacy first** — Never commit DNA data
2. **Not medical advice** — Always include disclaimers
3. **Single HTML output** — No external dependencies except Google Fonts
4. **Russian language** — Webpage output in Russian
5. **Highlight risks** — Use red for critical findings

