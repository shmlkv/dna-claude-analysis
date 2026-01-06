#!/usr/bin/env python3
"""
Pain Sensitivity SNP Analysis Script
Analyzes pain-related genetic markers from 23andMe data
"""

import os
from collections import defaultdict
from datetime import datetime

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Pain Sensitivity Categories
# =============================================================================

PAIN_SNPS = {
    "pain_threshold": {
        "name": "Базовая чувствительность к боли",
        "snps": {
            "rs4680": {
                "gene": "COMT Val158Met",
                "description": "Катехол-О-метилтрансфераза - метаболизм дофамина и катехоламинов",
                "risk_allele": "A",
                "interpretation": {
                    "GG": ("low_pain", "Val/Val - низкая чувствительность к боли, быстрый метаболизм катехоламинов"),
                    "AG": ("moderate", "Val/Met - средняя чувствительность к боли"),
                    "AA": ("high_pain", "Met/Met - высокая чувствительность к боли, медленный метаболизм катехоламинов"),
                }
            },
            "rs6746030": {
                "gene": "SCN9A",
                "description": "Натриевый канал Nav1.7 - проведение болевых сигналов",
                "risk_allele": "A",
                "interpretation": {
                    "GG": ("low_pain", "Более высокий болевой порог"),
                    "AG": ("moderate", "Средний болевой порог"),
                    "AA": ("high_pain", "Более низкий болевой порог, повышенная чувствительность"),
                }
            },
            "rs8007267": {
                "gene": "GCH1",
                "description": "ГТФ-циклогидролаза 1 - часть защитного гаплотипа",
                "risk_allele": "A",
                "interpretation": {
                    "GG": ("normal", "Стандартный вариант"),
                    "AG": ("protective", "Носитель защитного аллеля - сниженная болевая чувствительность"),
                    "AA": ("protective", "Защитный генотип - сниженная болевая чувствительность"),
                }
            },
            "rs3783641": {
                "gene": "GCH1",
                "description": "ГТФ-циклогидролаза 1 - защитный гаплотип против хронической боли",
                "risk_allele": "A",
                "interpretation": {
                    "TT": ("normal", "Стандартный вариант"),
                    "AT": ("protective", "Носитель защитного аллеля"),
                    "AA": ("protective", "Защитный генотип против хронической боли"),
                }
            },
            "rs10483639": {
                "gene": "GCH1",
                "description": "ГТФ-циклогидролаза 1 - компонент болезащитного гаплотипа",
                "risk_allele": "G",
                "interpretation": {
                    "CC": ("normal", "Стандартный вариант"),
                    "CG": ("protective", "Носитель защитного аллеля"),
                    "GG": ("protective", "Защитный генотип - меньше риск хронической боли"),
                }
            },
        }
    },

    "opioid_response": {
        "name": "Ответ на опиоидные анальгетики",
        "snps": {
            "rs1799971": {
                "gene": "OPRM1 A118G",
                "description": "Мю-опиоидный рецептор - основная мишень опиоидов",
                "risk_allele": "G",
                "interpretation": {
                    "AA": ("normal", "Нормальный ответ на опиоиды, стандартные дозы"),
                    "AG": ("moderate", "Сниженный ответ - могут потребоваться повышенные дозы"),
                    "GG": ("high", "Значительно сниженный ответ - нужны существенно выше дозы опиоидов"),
                }
            },
            "rs1045642": {
                "gene": "ABCB1 C3435T",
                "description": "P-гликопротеин - транспорт опиоидов через ГЭБ",
                "risk_allele": "T",
                "interpretation": {
                    "CC": ("good", "Хороший ответ на опиоиды, лучшее проникновение в ЦНС"),
                    "CT": ("moderate", "Средний ответ на опиоиды"),
                    "TT": ("poor", "Сниженный ответ на опиоиды, хуже проникновение"),
                }
            },
        }
    },

    "inflammatory_pain": {
        "name": "Воспалительная боль",
        "snps": {
            "rs1800629": {
                "gene": "TNF-α -308G>A",
                "description": "Фактор некроза опухоли альфа - ключевой медиатор воспаления",
                "risk_allele": "A",
                "interpretation": {
                    "GG": ("normal", "Нормальная продукция TNF-α"),
                    "AG": ("moderate", "Повышенная продукция TNF-α, склонность к воспалительной боли"),
                    "AA": ("high", "Высокая продукция TNF-α, повышен риск хронического воспаления и боли"),
                }
            },
            "rs1800795": {
                "gene": "IL-6 -174G>C",
                "description": "Интерлейкин-6 - провоспалительный цитокин",
                "risk_allele": "C",
                "interpretation": {
                    "GG": ("normal", "Нормальный уровень IL-6"),
                    "GC": ("moderate", "Умеренно повышенный IL-6"),
                    "CG": ("moderate", "Умеренно повышенный IL-6"),
                    "CC": ("high", "Повышенная продукция IL-6, склонность к воспалительной боли"),
                }
            },
            "rs1800896": {
                "gene": "IL-10 -1082A>G",
                "description": "Интерлейкин-10 - противовоспалительный цитокин",
                "risk_allele": "A",
                "interpretation": {
                    "GG": ("protective", "Высокая продукция IL-10, хорошая противовоспалительная защита"),
                    "AG": ("normal", "Средняя продукция IL-10"),
                    "GA": ("normal", "Средняя продукция IL-10"),
                    "AA": ("high", "Низкая продукция IL-10, сниженная противовоспалительная защита"),
                }
            },
        }
    },

    "migraine": {
        "name": "Мигрень",
        "snps": {
            "rs1835740": {
                "gene": "MTDH/AEG-1",
                "description": "Метадгерин - регуляция глутамата, первый GWAS-маркер мигрени",
                "risk_allele": "C",
                "interpretation": {
                    "TT": ("normal", "Нормальный риск мигрени"),
                    "CT": ("moderate", "Умеренно повышен риск мигрени (~20%)"),
                    "TC": ("moderate", "Умеренно повышен риск мигрени (~20%)"),
                    "CC": ("high", "Повышен риск мигрени (~40%)"),
                }
            },
            "rs2651899": {
                "gene": "PRDM16",
                "description": "Транскрипционный фактор - развитие нервной системы",
                "risk_allele": "C",
                "interpretation": {
                    "TT": ("normal", "Нормальный риск мигрени"),
                    "CT": ("moderate", "Немного повышен риск мигрени"),
                    "TC": ("moderate", "Немного повышен риск мигрени"),
                    "CC": ("high", "Повышен риск мигрени"),
                }
            },
            "rs10166942": {
                "gene": "TRPM8",
                "description": "Холодовой рецептор - сенсорная чувствительность",
                "risk_allele": "T",
                "interpretation": {
                    "CC": ("normal", "Нормальный риск мигрени"),
                    "CT": ("moderate", "Немного повышен риск мигрени"),
                    "TC": ("moderate", "Немного повышен риск мигрени"),
                    "TT": ("high", "Повышен риск мигрени с аурой"),
                }
            },
            "rs11172113": {
                "gene": "LRP1",
                "description": "Рецептор липопротеинов низкой плотности - нейропротекция",
                "risk_allele": "T",
                "interpretation": {
                    "CC": ("normal", "Нормальный риск мигрени"),
                    "CT": ("moderate", "Немного повышен риск"),
                    "TC": ("moderate", "Немного повышен риск"),
                    "TT": ("high", "Повышен риск мигрени"),
                }
            },
        }
    },

    "anesthetics": {
        "name": "Местные анестетики",
        "snps": {
            "rs12532": {
                "gene": "SCN5A",
                "description": "Натриевый канал - мишень местных анестетиков",
                "risk_allele": "T",
                "interpretation": {
                    "CC": ("normal", "Нормальный ответ на местные анестетики"),
                    "CT": ("moderate", "Немного сниженная эффективность анестезии"),
                    "TC": ("moderate", "Немного сниженная эффективность анестезии"),
                    "TT": ("high", "Может потребоваться больше анестетика"),
                }
            },
            "rs1805007": {
                "gene": "MC1R R151C",
                "description": "Меланокортиновый рецептор 1 - ассоциирован с рыжими волосами",
                "risk_allele": "T",
                "interpretation": {
                    "CC": ("normal", "Нормальный ответ на анестезию"),
                    "CT": ("high", "Носитель варианта MC1R - может требоваться на 19% больше анестетика"),
                    "TC": ("high", "Носитель варианта MC1R - может требоваться на 19% больше анестетика"),
                    "TT": ("high", "Вариант MC1R - требуется значительно больше местного анестетика"),
                }
            },
        }
    },

    "nsaids": {
        "name": "Ответ на НПВС (нестероидные противовоспалительные)",
        "snps": {
            "rs1799853": {
                "gene": "CYP2C9*2",
                "description": "Цитохром P450 2C9 - метаболизм НПВС",
                "risk_allele": "T",
                "interpretation": {
                    "CC": ("normal", "Нормальный метаболизатор - стандартные дозы НПВС"),
                    "CT": ("moderate", "Промежуточный метаболизатор - осторожность с дозой"),
                    "TC": ("moderate", "Промежуточный метаболизатор - осторожность с дозой"),
                    "TT": ("high", "Медленный метаболизатор - повышен риск побочных эффектов НПВС, снизить дозу"),
                }
            },
            "rs1057910": {
                "gene": "CYP2C9*3",
                "description": "Цитохром P450 2C9 - важнейший фермент метаболизма НПВС",
                "risk_allele": "C",
                "interpretation": {
                    "AA": ("normal", "Нормальный метаболизатор - стандартные дозы"),
                    "AC": ("moderate", "Промежуточный метаболизатор - риск накопления НПВС"),
                    "CA": ("moderate", "Промежуточный метаболизатор - риск накопления НПВС"),
                    "CC": ("high", "Медленный метаболизатор - высокий риск побочных эффектов, снизить дозу на 50%"),
                }
            },
        }
    },
}


def load_genome():
    """Load genome data into a dictionary"""
    genome = {}
    with open(GENOME_FILE, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
                genome[rsid] = {
                    'chromosome': chrom,
                    'position': pos,
                    'genotype': genotype
                }
    return genome


def normalize_genotype(genotype):
    """Normalize genotype for comparison (sort alleles)"""
    if len(genotype) == 2:
        return ''.join(sorted(genotype))
    return genotype


def analyze_snp(snp_id, snp_info, genome_data):
    """Analyze a single SNP"""
    result = {
        'snp_id': snp_id,
        'gene': snp_info['gene'],
        'description': snp_info['description'],
        'risk_allele': snp_info['risk_allele'],
        'found': False,
        'genotype': None,
        'risk_level': None,
        'interpretation': None
    }

    if snp_id in genome_data:
        result['found'] = True
        raw_genotype = genome_data[snp_id]['genotype']
        result['genotype'] = raw_genotype
        result['chromosome'] = genome_data[snp_id]['chromosome']
        result['position'] = genome_data[snp_id]['position']

        # Try to find interpretation
        normalized = normalize_genotype(raw_genotype)
        interpretations = snp_info.get('interpretation', {})

        # Try both original and normalized genotype
        for gt in [raw_genotype, normalized]:
            if gt in interpretations:
                result['risk_level'], result['interpretation'] = interpretations[gt]
                break

        # If still not found, try reverse
        if result['interpretation'] is None and len(raw_genotype) == 2:
            reversed_gt = raw_genotype[::-1]
            if reversed_gt in interpretations:
                result['risk_level'], result['interpretation'] = interpretations[reversed_gt]

    return result


def analyze_gch1_haplotype(genome):
    """Analyze GCH1 protective haplotype"""
    rs8007267 = genome.get('rs8007267', {}).get('genotype', '')
    rs3783641 = genome.get('rs3783641', {}).get('genotype', '')
    rs10483639 = genome.get('rs10483639', {}).get('genotype', '')

    protective_count = 0
    details = []

    if 'A' in rs8007267:
        protective_count += rs8007267.count('A')
        details.append(f"rs8007267: {rs8007267}")
    if 'A' in rs3783641:
        protective_count += rs3783641.count('A')
        details.append(f"rs3783641: {rs3783641}")
    if 'G' in rs10483639:
        protective_count += rs10483639.count('G')
        details.append(f"rs10483639: {rs10483639}")

    if protective_count >= 4:
        status = ('strong_protective', 'Сильный защитный GCH1 гаплотип - значительно снижена болевая чувствительность')
    elif protective_count >= 2:
        status = ('protective', 'Частичный защитный GCH1 гаплотип - умеренно снижена болевая чувствительность')
    elif protective_count >= 1:
        status = ('mild_protective', 'Носитель защитных аллелей GCH1')
    else:
        status = ('normal', 'Нет защитного GCH1 гаплотипа')

    return {
        'rs8007267': rs8007267,
        'rs3783641': rs3783641,
        'rs10483639': rs10483639,
        'protective_alleles': protective_count,
        'status': status[0],
        'interpretation': status[1],
        'details': details
    }


def analyze_cyp2c9_status(results):
    """Determine combined CYP2C9 metabolizer status for NSAIDs"""
    cyp2c9_2 = None
    cyp2c9_3 = None

    for r in results:
        if r['snp_id'] == 'rs1799853':
            cyp2c9_2 = r['genotype']
        elif r['snp_id'] == 'rs1057910':
            cyp2c9_3 = r['genotype']

    if not cyp2c9_2 and not cyp2c9_3:
        return None

    # Count variant alleles
    variant_count = 0
    if cyp2c9_2:
        variant_count += cyp2c9_2.count('T')
    if cyp2c9_3:
        variant_count += cyp2c9_3.count('C')

    if variant_count >= 3:
        status = ('poor', 'Плохой метаболизатор CYP2C9 - снизить дозу НПВС на 50-75%')
    elif variant_count == 2:
        status = ('intermediate', 'Промежуточный метаболизатор CYP2C9 - снизить дозу НПВС на 25-50%')
    elif variant_count == 1:
        status = ('intermediate', 'Промежуточный метаболизатор CYP2C9 - осторожность с дозой')
    else:
        status = ('normal', 'Нормальный метаболизатор CYP2C9 - стандартные дозы НПВС')

    return {
        'cyp2c9_2': cyp2c9_2,
        'cyp2c9_3': cyp2c9_3,
        'variant_alleles': variant_count,
        'status': status[0],
        'interpretation': status[1]
    }


def analyze_comt_pain_profile(genome):
    """Analyze COMT Val158Met for pain sensitivity profile"""
    rs4680 = genome.get('rs4680', {}).get('genotype', '')

    if rs4680 == 'GG':
        profile = {
            'type': 'Warrior (Val/Val)',
            'pain_sensitivity': 'Низкая',
            'stress_tolerance': 'Высокая',
            'dopamine_metabolism': 'Быстрый',
            'clinical': 'Меньше нужны обезболивающие, лучше переносит стресс'
        }
    elif rs4680 == 'AA':
        profile = {
            'type': 'Worrier (Met/Met)',
            'pain_sensitivity': 'Высокая',
            'stress_tolerance': 'Низкая',
            'dopamine_metabolism': 'Медленный',
            'clinical': 'Может потребоваться больше обезболивающих, выше когнитивные функции'
        }
    elif rs4680 in ['AG', 'GA']:
        profile = {
            'type': 'Смешанный (Val/Met)',
            'pain_sensitivity': 'Средняя',
            'stress_tolerance': 'Средняя',
            'dopamine_metabolism': 'Средний',
            'clinical': 'Сбалансированный профиль боли и стресса'
        }
    else:
        profile = {
            'type': 'Не определён',
            'pain_sensitivity': 'Н/Д',
            'stress_tolerance': 'Н/Д',
            'dopamine_metabolism': 'Н/Д',
            'clinical': 'Генотип не найден'
        }

    profile['genotype'] = rs4680
    return profile


def generate_category_report(category, results, genome):
    """Generate report for a category"""
    cat_info = PAIN_SNPS[category]

    report = []
    report.append(f"# {cat_info['name']}")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("\n## Результаты\n")

    # Statistics
    found = sum(1 for r in results if r['found'])
    report.append(f"Найдено маркеров: {found}/{len(results)}\n")

    # Risk summary
    risk_counts = defaultdict(int)
    for r in results:
        if r['risk_level']:
            risk_counts[r['risk_level']] += 1

    if risk_counts:
        report.append("### Сводка по результатам\n")
        risk_emoji = {
            'high': '🔴',
            'high_pain': '🔴',
            'moderate': '🟡',
            'low': '🟢',
            'low_pain': '🟢',
            'normal': '✅',
            'protective': '🛡️',
            'strong_protective': '🛡️🛡️',
            'mild_protective': '🛡️',
            'good': '✅',
            'poor': '🔴',
            'info': 'ℹ️'
        }
        for risk, count in sorted(risk_counts.items()):
            emoji = risk_emoji.get(risk, '•')
            report.append(f"- {emoji} {risk}: {count}")

    report.append("\n### Детальные результаты\n")
    report.append("| SNP | Ген | Генотип | Статус | Интерпретация |")
    report.append("|-----|-----|---------|--------|---------------|")

    for r in results:
        if r['found']:
            risk_label = r['risk_level'] or 'н/д'
            interp = r['interpretation'] or 'Нет данных'
            report.append(f"| {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {risk_label} | {interp} |")
        else:
            report.append(f"| {r['snp_id']} | {r['gene']} | - | - | Не найден в геноме |")

    # Special sections for specific categories
    if category == 'pain_threshold':
        report.append("\n### GCH1 защитный гаплотип\n")
        gch1 = analyze_gch1_haplotype(genome)
        report.append(f"- rs8007267: {gch1['rs8007267'] or 'не найден'}")
        report.append(f"- rs3783641: {gch1['rs3783641'] or 'не найден'}")
        report.append(f"- rs10483639: {gch1['rs10483639'] or 'не найден'}")
        report.append(f"- Защитных аллелей: {gch1['protective_alleles']}")
        report.append(f"- **Статус: {gch1['status']}**")
        report.append(f"- {gch1['interpretation']}")

        report.append("\n### COMT профиль боли\n")
        comt = analyze_comt_pain_profile(genome)
        report.append(f"- Генотип rs4680: **{comt['genotype']}**")
        report.append(f"- Тип: **{comt['type']}**")
        report.append(f"- Болевая чувствительность: {comt['pain_sensitivity']}")
        report.append(f"- Устойчивость к стрессу: {comt['stress_tolerance']}")
        report.append(f"- Метаболизм дофамина: {comt['dopamine_metabolism']}")
        report.append(f"- Клиническое значение: {comt['clinical']}")

    if category == 'nsaids':
        cyp2c9 = analyze_cyp2c9_status(results)
        if cyp2c9:
            report.append("\n### Статус CYP2C9 метаболизатора\n")
            report.append(f"- CYP2C9*2 (rs1799853): {cyp2c9['cyp2c9_2'] or 'не найден'}")
            report.append(f"- CYP2C9*3 (rs1057910): {cyp2c9['cyp2c9_3'] or 'не найден'}")
            report.append(f"- Вариантных аллелей: {cyp2c9['variant_alleles']}")
            report.append(f"- **Статус: {cyp2c9['status']}**")
            report.append(f"- {cyp2c9['interpretation']}")

    return '\n'.join(report)


def generate_summary_report(all_results, genome):
    """Generate overall summary report"""
    report = []
    report.append("# Анализ чувствительности к боли")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("\n---\n")

    report.append("## Важные предупреждения\n")
    report.append("1. **Это НЕ медицинский диагноз** - только информационный генетический анализ")
    report.append("2. **Болевая чувствительность зависит от многих факторов** - генетика лишь часть")
    report.append("3. **Индивидуальный ответ может отличаться** - генотип не гарантирует фенотип")
    report.append("4. **Для назначения лекарств** - консультация врача обязательна\n")

    report.append("---\n")

    # COMT Pain Profile - Main highlight
    report.append("## Основной профиль боли (COMT)\n")
    comt = analyze_comt_pain_profile(genome)
    report.append(f"| Параметр | Значение |")
    report.append(f"|----------|----------|")
    report.append(f"| Генотип | **{comt['genotype']}** |")
    report.append(f"| Тип | **{comt['type']}** |")
    report.append(f"| Болевая чувствительность | {comt['pain_sensitivity']} |")
    report.append(f"| Устойчивость к стрессу | {comt['stress_tolerance']} |")
    report.append(f"| Метаболизм дофамина | {comt['dopamine_metabolism']} |")
    report.append(f"\n**Клиническое значение:** {comt['clinical']}\n")

    # GCH1 Haplotype
    report.append("## Защитный гаплотип GCH1\n")
    gch1 = analyze_gch1_haplotype(genome)
    if gch1['protective_alleles'] > 0:
        report.append(f"- Статус: **{gch1['status']}**")
        report.append(f"- {gch1['interpretation']}")
        report.append(f"- Защитных аллелей: {gch1['protective_alleles']}/6")
    else:
        report.append("- Защитный гаплотип не обнаружен")
    report.append("")

    report.append("---\n")

    # Collect findings by importance
    high_pain = []
    low_pain = []
    medication_warnings = []
    protective = []

    for category, results in all_results.items():
        cat_name = PAIN_SNPS[category]['name']
        for r in results:
            if r['risk_level'] in ['high', 'high_pain', 'poor']:
                if category in ['opioid_response', 'nsaids', 'anesthetics']:
                    medication_warnings.append((cat_name, r))
                else:
                    high_pain.append((cat_name, r))
            elif r['risk_level'] in ['low_pain', 'good']:
                low_pain.append((cat_name, r))
            elif r['risk_level'] in ['protective', 'strong_protective', 'mild_protective']:
                protective.append((cat_name, r))

    if medication_warnings:
        report.append("## Важно для приёма лекарств\n")
        report.append("| Категория | SNP | Ген | Генотип | Значение |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in medication_warnings:
            report.append(f"| {cat} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if high_pain:
        report.append("## Повышенная болевая чувствительность\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in high_pain:
            report.append(f"| {cat} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if low_pain:
        report.append("## Сниженная болевая чувствительность\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in low_pain:
            report.append(f"| {cat} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if protective:
        report.append("## Защитные варианты\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in protective:
            report.append(f"| {cat} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    # CYP2C9 Status for NSAIDs
    nsaid_results = all_results.get('nsaids', [])
    cyp2c9 = analyze_cyp2c9_status(nsaid_results)
    if cyp2c9:
        report.append("## Метаболизм НПВС (CYP2C9)\n")
        report.append(f"- CYP2C9*2: {cyp2c9['cyp2c9_2'] or 'не найден'}")
        report.append(f"- CYP2C9*3: {cyp2c9['cyp2c9_3'] or 'не найден'}")
        report.append(f"- **Статус метаболизатора: {cyp2c9['status']}**")
        report.append(f"- Рекомендация: {cyp2c9['interpretation']}\n")

    # Migraine risk
    report.append("## Риск мигрени\n")
    migraine_results = all_results.get('migraine', [])
    migraine_risk_count = sum(1 for r in migraine_results if r['risk_level'] in ['high', 'moderate'])
    if migraine_risk_count >= 3:
        report.append("- **Повышенный генетический риск мигрени**")
    elif migraine_risk_count >= 1:
        report.append("- Умеренный генетический риск мигрени")
    else:
        report.append("- Нет повышенного генетического риска мигрени")

    migraine_found = [r for r in migraine_results if r['found']]
    if migraine_found:
        report.append("\nРелевантные маркеры:")
        for r in migraine_found:
            status = r['risk_level'] or 'н/д'
            report.append(f"- {r['snp_id']} ({r['gene']}): {r['genotype']} - {status}")
    report.append("")

    report.append("---\n")
    report.append("## Статистика анализа\n")

    total_snps = 0
    found_snps = 0
    for cat, results in all_results.items():
        total_snps += len(results)
        found_snps += sum(1 for r in results if r['found'])

    report.append(f"- Всего проанализировано SNP: {total_snps}")
    report.append(f"- Найдено в геноме: {found_snps}")
    report.append(f"- Не найдено: {total_snps - found_snps}")

    report.append("\n---\n")
    report.append("## Практические рекомендации\n")

    # Generate personalized recommendations
    recommendations = []

    if comt['pain_sensitivity'] == 'Высокая':
        recommendations.append("- **Болевая чувствительность:** Вы можете быть более чувствительны к боли. Обсудите с врачом индивидуальный подход к обезболиванию.")

    if comt['pain_sensitivity'] == 'Низкая':
        recommendations.append("- **Болевая чувствительность:** У вас, вероятно, высокий болевой порог. Помните, что боль - важный сигнал организма.")

    opioid_results = all_results.get('opioid_response', [])
    for r in opioid_results:
        if r['snp_id'] == 'rs1799971' and r['risk_level'] in ['moderate', 'high']:
            recommendations.append("- **Опиоиды:** Возможно сниженный ответ на опиоидные анальгетики. Сообщите анестезиологу.")
            break

    if cyp2c9 and cyp2c9['status'] != 'normal':
        recommendations.append(f"- **НПВС:** {cyp2c9['interpretation']}")

    anesthetic_results = all_results.get('anesthetics', [])
    for r in anesthetic_results:
        if r['snp_id'] == 'rs1805007' and r['risk_level'] == 'high':
            recommendations.append("- **Местная анестезия:** Может потребоваться повышенная доза местных анестетиков. Предупредите стоматолога/хирурга.")
            break

    if migraine_risk_count >= 2:
        recommendations.append("- **Мигрень:** Повышенная генетическая предрасположенность. Избегайте известных триггеров, ведите дневник головной боли.")

    if recommendations:
        report.extend(recommendations)
    else:
        report.append("- Специфических рекомендаций по результатам анализа нет.")

    return '\n'.join(report)


def main():
    print("=" * 60)
    print("АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ К БОЛИ ПО ГЕНОМУ")
    print("=" * 60)

    print("\n[1/4] Загрузка генома...")
    genome = load_genome()
    print(f"      Загружено {len(genome)} SNP")

    print("\n[2/4] Анализ маркеров по категориям...")
    all_results = {}

    for category, cat_info in PAIN_SNPS.items():
        print(f"      -> {cat_info['name']}...")
        results = []
        for snp_id, snp_info in cat_info['snps'].items():
            result = analyze_snp(snp_id, snp_info, genome)
            results.append(result)
        all_results[category] = results

        # Count found
        found = sum(1 for r in results if r['found'])
        print(f"        Найдено: {found}/{len(results)}")

    print("\n[3/4] Генерация детальных отчётов...")
    for category, results in all_results.items():
        report = generate_category_report(category, results, genome)
        report_dir = f"{REPORTS_PATH}/pain/{category}"
        os.makedirs(report_dir, exist_ok=True)
        report_path = f"{report_dir}/report.md"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"      -> {report_path}")

    print("\n[4/4] Генерация сводного отчёта...")
    summary = generate_summary_report(all_results, genome)
    summary_path = f"{REPORTS_PATH}/pain/report.md"
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(summary)
    print(f"      -> {summary_path}")

    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)

    # Print key findings to console
    print("\nКЛЮЧЕВЫЕ НАХОДКИ:\n")

    # COMT Profile
    comt = analyze_comt_pain_profile(genome)
    print(f"COMT профиль: {comt['type']}")
    print(f"  - Болевая чувствительность: {comt['pain_sensitivity']}")
    print(f"  - {comt['clinical']}")
    print()

    # Medication warnings
    for category, results in all_results.items():
        warnings = [r for r in results if r['risk_level'] in ['high', 'poor']]
        if warnings and category in ['opioid_response', 'nsaids', 'anesthetics']:
            print(f"ВНИМАНИЕ - {PAIN_SNPS[category]['name']}:")
            for r in warnings:
                print(f"    * {r['gene']} ({r['genotype']}): {r['interpretation']}")
            print()


if __name__ == "__main__":
    main()
