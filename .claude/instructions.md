# Webpage Generation Instructions

## When to Generate

User asks for:
- "Generate webpage"
- "Create DNA Terminal"
- "Build HTML from reports"
- "Визуализируй данные"

## Generation Steps

### 1. Read Reports
```
reports/
├── health_report.md
├── ancestry_report.md
├── nutrition_report.md
└── ... (all *_report.md files)
```

### 2. Parse Data
Extract from each report:
- Gene/SNP names
- Genotypes (AA, AG, GG, etc.)
- Status (normal, risk, moderate)
- Descriptions

### 3. Build HTML
Follow `webpage/STYLE_GUIDE.md`, look at `example.html` for reference:
- Single file with inline styles
- JetBrains Mono font
- Terminal color palette
- Fixed navigation header

### 4. Create Sections
Map reports to sections:
| Report | Section ID | Nav Label |
|--------|------------|-----------|
| health_report.md | health | HEALTH |
| ancestry_report.md | ancestry | ANCESTRY |
| nutrition_report.md | nutrition | NUTRITION |
| sports_fitness_report.md | sports | SPORTS |
| ... | ... | ... |

### 5. Color Coding
Apply status colors:
- `status-normal` → green (#00ff41)
- `status-moderate` → amber (#ffb000)
- `status-risk` → red (#ff3333)
- `status-info` → cyan (#00ffff)
- `status-na` → gray (#666666)

### 6. Output
Save to: `webpage/dna_terminal.html`

## Content Language
Russian — all text content in Russian.

## Required Elements

### Header
```html
<header>
    <div class="header-top">
        <div class="logo">DNA TERMINAL v1.0</div>
        <div class="timestamp">[DATE] | [USERNAME]</div>
    </div>
    <nav>[SECTION LINKS]</nav>
</header>
```

### Status Box
First section with:
- Genome info (source, SNP count)
- Key alerts (risks in red)
- Key highlights (advantages in green)
- Analysis progress bar

### Disclaimers
Include in health-related sections:
```html
<div class="disclaimer">
    <strong>! Важно</strong>
    Данные носят информационный характер...
</div>
```

## Testing
After generation, offer to open in browser for review.

