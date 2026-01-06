# DNA TERMINAL — Style & Structure Guide

## Overview

Single-page web application styled as a **terminal/console interface** for genetic data visualization. Uses hacker/cyberpunk aesthetic with characteristic green text on dark background.

---

## Document Structure

```
<!DOCTYPE html>
├── <html lang="ru">
│   ├── <head>
│   │   ├── meta (charset, viewport)
│   │   ├── title
│   │   ├── Google Fonts (JetBrains Mono)
│   │   └── <style> (all CSS inline)
│   │
│   └── <body>
│       ├── <header> (fixed header)
│       │   ├── .header-top (logo + timestamp)
│       │   └── <nav> (anchor links)
│       │
│       ├── <main>
│       │   ├── .status-box (general info)
│       │   └── <section> × N (data sections)
│       │
│       └── <script> (smooth scroll + active nav)
```

---

## Color Palette

```css
:root {
    --bg: #0a0a0a;           /* Main background */
    --bg-light: #111111;      /* Element background */
    --bg-section: #0d0d0d;    /* Section background */
    --green: #00ff41;         /* Primary accent (Matrix green) */
    --green-dim: #00cc33;     /* Dimmed green */
    --amber: #ffb000;         /* Warnings */
    --red: #ff3333;           /* Risks/errors */
    --cyan: #00ffff;          /* Info accent */
    --gray: #666666;          /* Borders, dividers */
    --gray-light: #888888;    /* Secondary text */
    --white: #e0e0e0;         /* Primary data text */
}
```

### Color Usage

| Color | Usage |
|-------|-------|
| `--green` | Headers, accents, positive statuses |
| `--amber` | Moderate risks, warnings |
| `--red` | High risks, critical data |
| `--cyan` | Info, data highlights, genotype values |
| `--gray-light` | Labels, secondary text |
| `--white` | Table values |

---

## Typography

### Font
```html
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;700&display=swap" rel="stylesheet">
```

```css
body {
    font-family: 'JetBrains Mono', monospace;
    font-size: 14px;
    line-height: 1.6;
}
```

### Size Hierarchy

| Element | Size | Weight | Extra |
|---------|------|--------|-------|
| Logo | 18px | 700 | — |
| Section header | 14px | — | uppercase, letter-spacing: 2px |
| Subsection title | 12px | — | — |
| Navigation | 11px | — | uppercase, letter-spacing: 1px |
| Tables | 12px | — | — |
| Table headers | 10px | 500 | uppercase, letter-spacing: 1px |
| Score value | 48px | 700 | — |
| Profile value | 18px | 700 | — |

---

## Components

### 1. Header (fixed)

```html
<header>
    <div class="header-top">
        <div class="logo">DNA TERMINAL v1.0</div>
        <div class="timestamp">2026-01-06 16:29 | Username</div>
    </div>
    <nav>
        <a href="#section-id">SECTION_NAME</a>
        <!-- ... -->
    </nav>
</header>
```

**Features:**
- `position: fixed` with `z-index: 1000`
- Bottom border: `1px solid var(--green-dim)`
- Logo with animated cursor:
```css
.logo::before {
    content: ">";
    margin-right: 8px;
    animation: blink 1s step-end infinite;
}
```

### 2. Status Box (main block)

```html
<div class="status-box" id="status">
    <div class="status-header">// SYSTEM STATUS</div>
    <div class="status-grid">
        <div class="status-item">
            <span class="status-label">LABEL:</span>
            <span class="status-value">Value</span>
        </div>
    </div>
    <div class="alerts-grid">
        <div class="alerts-column warnings">
            <h3>// ALERTS</h3>
            <div class="alert-item warning">
                <span class="alert-icon">!</span>
                <span>Alert text</span>
            </div>
        </div>
        <div class="alerts-column highlights">
            <h3>// HIGHLIGHTS</h3>
            <div class="alert-item success">
                <span class="alert-icon">+</span>
                <span>Highlight text</span>
            </div>
        </div>
    </div>
    <div class="progress-bar">
        <div class="progress-track">
            <div class="progress-fill" style="flex: 37;"></div>
            <div class="progress-empty" style="flex: 4;"></div>
        </div>
        <div class="progress-label">37/41 items analyzed</div>
    </div>
</div>
```

### 3. Section (data section)

```html
<section id="section-id">
    <div class="section-header">
        <span class="section-title">// SECTION NAME</span>
        <span class="section-badge good|moderate|risk">BADGE TEXT</span>
    </div>
    <div class="section-content">
        <!-- section content -->
    </div>
</section>
```

**Status badges:**
```css
.section-badge.good { background: rgba(0, 255, 65, 0.2); color: var(--green); }
.section-badge.moderate { background: rgba(255, 176, 0, 0.2); color: var(--amber); }
.section-badge.risk { background: rgba(255, 51, 51, 0.2); color: var(--red); }
```

### 4. Subsection

```html
<div class="subsection">
    <div class="subsection-title">SUBSECTION NAME</div>
    <!-- table or other content -->
</div>
```

### 5. Table

```html
<table>
    <tr>
        <th>Column 1</th>
        <th>Column 2</th>
        <th>Status</th>
    </tr>
    <tr>
        <td>Value</td>
        <td class="geno">AA</td>
        <td class="status-normal|moderate|risk|info|na">STATUS TEXT</td>
    </tr>
</table>
```

**Status classes:**
```css
.geno { font-weight: 700; color: var(--cyan); }
.status-normal { color: var(--green); }
.status-moderate { color: var(--amber); }
.status-risk { color: var(--red); }
.status-info { color: var(--cyan); }
.status-na { color: var(--gray); }
```

### 6. Bar Chart

```html
<div class="bar-container">
    <div class="bar-label">
        <span>Label</span>
        <span class="status-normal">82%</span>
    </div>
    <div class="bar-track">
        <div class="bar-fill" style="width: 82%;"></div>
    </div>
</div>
```

**Color modifiers:**
```css
.bar-fill { background: var(--green); }
.bar-fill.amber { background: var(--amber); }
.bar-fill.cyan { background: var(--cyan); }
.bar-fill.red { background: var(--red); }
```

### 7. Finding (finding card)

```html
<div class="finding warning|risk|good|info">
    <div class="finding-icon">!</div>
    <div class="finding-content">
        <div class="finding-gene">GENE_NAME (Genotype)</div>
        <div class="finding-desc">Description text</div>
    </div>
</div>
```

**Border colors:**
```css
.finding { border-left: 3px solid var(--gray); }
.finding.warning { border-left-color: var(--amber); }
.finding.risk { border-left-color: var(--red); }
.finding.good { border-left-color: var(--green); }
.finding.info { border-left-color: var(--cyan); }
```

### 8. Profile Box

```html
<div class="profile-box">
    <div class="profile-title">Title</div>
    <div class="profile-value" style="color: var(--cyan);">Value</div>
    <div class="profile-desc">Description</div>
</div>
```

### 9. Score Display

```html
<div class="score-display">
    <div class="score-value">97</div>
    <div class="score-meta">
        <strong>Percentile</strong>
        Additional info
    </div>
</div>
```

### 10. Disclaimer

```html
<div class="disclaimer">
    <strong>! Important</strong>
    Warning text content.
</div>
```

**Critical variant:**
```html
<div class="disclaimer" style="background: rgba(255, 51, 51, 0.1); border-color: var(--red); color: var(--red);">
    <strong>!! CRITICAL</strong>
    Critical warning content.
</div>
```

### 11. Two Columns Layout

```html
<div class="two-cols">
    <div>Column 1 content</div>
    <div>Column 2 content</div>
</div>
```

---

## Responsiveness

```css
@media (max-width: 768px) {
    .alerts-grid { grid-template-columns: 1fr; }
    .two-cols { grid-template-columns: 1fr; }
}
```

---

## JavaScript Features

### Smooth Scroll

```javascript
document.querySelectorAll('nav a').forEach(anchor => {
    anchor.addEventListener('click', function(e) {
        e.preventDefault();
        const targetId = this.getAttribute('href');
        const target = document.querySelector(targetId);
        if (target) {
            const headerHeight = document.querySelector('header').offsetHeight;
            const targetPosition = target.getBoundingClientRect().top + window.pageYOffset - headerHeight - 20;
            window.scrollTo({
                top: targetPosition,
                behavior: 'smooth'
            });
        }
    });
});
```

### Active Navigation Highlight

```javascript
const sections = document.querySelectorAll('section, .status-box');
const navLinks = document.querySelectorAll('nav a');

function updateActiveNav() {
    const headerHeight = document.querySelector('header').offsetHeight;
    let current = '';

    sections.forEach(section => {
        const sectionTop = section.offsetTop - headerHeight - 100;
        if (window.pageYOffset >= sectionTop) {
            current = section.getAttribute('id');
        }
    });

    navLinks.forEach(link => {
        link.classList.remove('active');
        if (link.getAttribute('href') === '#' + current) {
            link.classList.add('active');
        }
    });
}

window.addEventListener('scroll', updateActiveNav);
updateActiveNav();
```

---

## Custom Scrollbar

```css
::-webkit-scrollbar { width: 8px; height: 8px; }
::-webkit-scrollbar-track { background: var(--bg); }
::-webkit-scrollbar-thumb { background: var(--gray); border-radius: 0; }
::-webkit-scrollbar-thumb:hover { background: var(--green-dim); }
```

---

## Step-by-Step Build Instructions

### Step 1: Base HTML

```html
<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DNA TERMINAL v1.0</title>
    <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;700&display=swap" rel="stylesheet">
    <style>
        /* CSS variables and styles */
    </style>
</head>
<body>
    <!-- content -->
    <script>
        /* JavaScript */
    </script>
</body>
</html>
```

### Step 2: CSS Variables

Add all CSS variables in `:root` for consistent colors.

### Step 3: Reset and Base Styles

```css
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
    font-family: 'JetBrains Mono', monospace;
    background: var(--bg);
    color: var(--green);
    font-size: 14px;
    line-height: 1.6;
}
```

### Step 4: Create Header

1. Add fixed header
2. Create logo and timestamp block
3. Add navigation with anchor links

### Step 5: Create Main

```html
<main>
    <!-- padding-top should account for header height -->
</main>
```

### Step 6: Add Status Box

First element — general genome information.

### Step 7: Add Sections

Each section:
1. `<section id="unique-id">`
2. `.section-header` with title and badge
3. `.section-content` with data

### Step 8: Populate Sections

Use component combinations:
- Tables for structured data
- Profile boxes for key metrics
- Finding cards for important discoveries
- Bar charts for percentages
- Two-cols for two-column layout

### Step 9: Add Navigation to Header

For each section add a link:
```html
<a href="#section-id">SECTION NAME</a>
```

### Step 10: Add JavaScript

Connect smooth scroll and active navigation.

---

## Sections List (standard set)

| ID | Name | Description |
|----|------|-------------|
| status | STATUS | General info, alerts |
| ancestry | ANCESTRY | Origin, haplogroups |
| health | HEALTH | Disease risks |
| nutrition | NUTRITION | Metabolism, food intolerances |
| sports | SPORTS | Athletic potential |
| psychology | PSYCH | Psychological traits |
| cognitive | COGNITIVE | Cognitive abilities |
| longevity | LONGEVITY | Longevity markers |
| sleep | SLEEP | Sleep, circadian rhythms |
| immunity | IMMUNITY | Immune system |
| pain | PAIN | Pain sensitivity |
| detox | DETOX | Detoxification |
| skin | SKIN | Skin |
| vision | VISION | Vision and hearing |
| traits | TRAITS | Physical traits |
| pharma | PHARMA | Pharmacogenomics |
| carrier | CARRIER | Carrier status |

---

## Design Principles

1. **Console aesthetic** — everything should look like a terminal
2. **Information density** — lots of data, but structured
3. **Color coding** — green/yellow/red for statuses
4. **Monospace font** — JetBrains Mono only
5. **Borders and frames** — clear block separation
6. **Comment-style** — headers in `// NAME` format
7. **Uppercase** — headers and labels in capitals
8. **Minimal rounding** — sharp corners, terminal look

---

## Pre-publish Checklist

- [ ] All sections have unique `id`
- [ ] Navigation contains all sections
- [ ] CSS variables used everywhere for colors
- [ ] Mobile responsiveness works
- [ ] Smooth scroll functions
- [ ] Active nav highlights on scroll
- [ ] All tables have headers
- [ ] Disclaimer present in medical sections
- [ ] Critical findings highlighted in red
