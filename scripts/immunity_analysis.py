#!/usr/bin/env python3
"""
Immunity SNP Analysis Script
Analyzes immunity-related genetic markers from 23andMe data
"""

import os
from collections import defaultdict
from datetime import datetime

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Organized by immunity category
# =============================================================================

IMMUNITY_SNPS = {
    "hla_system": {
        "name": "HLA-система (HLA-B27, целиакия DQ2/DQ8)",
        "snps": {
            "rs4349859": {
                "gene": "HLA-B27",
                "description": "Анкилозирующий спондилит, реактивный артрит",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "HLA-B27 положительный - высокий риск анкилозирующего спондилита"),
                    "AG": ("moderate", "Носитель HLA-B27 - повышенный риск спондилоартропатий"),
                    "GG": ("normal", "HLA-B27 отрицательный - нормальный риск"),
                }
            },
            "rs2187668": {
                "gene": "HLA-DQ2.5",
                "description": "Целиакия (главный маркер)",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "HLA-DQ2.5 гомозигота - очень высокий риск целиакии"),
                    "CT": ("moderate", "Носитель HLA-DQ2.5 - повышенный риск целиакии"),
                    "TC": ("moderate", "Носитель HLA-DQ2.5 - повышенный риск целиакии"),
                    "CC": ("normal", "Низкий риск целиакии по DQ2.5"),
                }
            },
            "rs7454108": {
                "gene": "HLA-DQ8",
                "description": "Целиакия, диабет 1 типа",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "HLA-DQ8 положительный - риск целиакии и СД1"),
                    "CT": ("low", "Носитель HLA-DQ8"),
                    "TC": ("low", "Носитель HLA-DQ8"),
                    "TT": ("normal", "Низкий риск по HLA-DQ8"),
                }
            },
        }
    },

    "autoimmune": {
        "name": "Аутоиммунные риски",
        "snps": {
            "rs2476601": {
                "gene": "PTPN22",
                "description": "Ревматоидный артрит, СД1, СКВ, тиреоидит",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Гомозигота риска - очень высокий риск аутоиммунных заболеваний"),
                    "AG": ("moderate", "Гетерозигота - повышенный риск RA, СД1, СКВ, тиреоидита"),
                    "GA": ("moderate", "Гетерозигота - повышенный риск RA, СД1, СКВ, тиреоидита"),
                    "GG": ("normal", "Нормальный риск аутоиммунных заболеваний"),
                }
            },
            "rs6457617": {
                "gene": "HLA-DRB1",
                "description": "Ревматоидный артрит (shared epitope)",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Высокий риск ревматоидного артрита"),
                    "CT": ("moderate", "Повышенный риск ревматоидного артрита"),
                    "TC": ("moderate", "Повышенный риск ревматоидного артрита"),
                    "CC": ("normal", "Нормальный риск ревматоидного артрита"),
                }
            },
            "rs3135388": {
                "gene": "HLA-DRB1 (MS)",
                "description": "Рассеянный склероз",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Высокий риск рассеянного склероза"),
                    "AG": ("moderate", "Повышенный риск рассеянного склероза"),
                    "GA": ("moderate", "Повышенный риск рассеянного склероза"),
                    "GG": ("normal", "Нормальный риск рассеянного склероза"),
                }
            },
            "rs2066847": {
                "gene": "NOD2",
                "description": "Болезнь Крона",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Высокий риск болезни Крона"),
                    "CG": ("moderate", "Повышенный риск болезни Крона (2-4x)"),
                    "GC": ("moderate", "Повышенный риск болезни Крона (2-4x)"),
                    "--": ("normal", "Инсерция отсутствует - нормальный риск"),
                    "GG": ("normal", "Нормальный риск болезни Крона"),
                }
            },
            "rs3087243": {
                "gene": "CTLA4",
                "description": "Аутоиммунный тиреоидит, СД1, RA",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("moderate", "Повышенный риск аутоиммунных заболеваний"),
                    "AG": ("low", "Немного повышенный риск"),
                    "GA": ("low", "Немного повышенный риск"),
                    "AA": ("normal", "Нормальный риск"),
                }
            },
        }
    },

    "cytokines": {
        "name": "Цитокины (воспалительный ответ)",
        "snps": {
            "rs1800629": {
                "gene": "TNF-alpha",
                "description": "Фактор некроза опухоли, воспаление",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Высокая продукция TNF-alpha - склонность к хроническому воспалению"),
                    "AG": ("moderate", "Повышенная продукция TNF-alpha"),
                    "GA": ("moderate", "Повышенная продукция TNF-alpha"),
                    "GG": ("normal", "Нормальная продукция TNF-alpha"),
                }
            },
            "rs1800795": {
                "gene": "IL-6",
                "description": "Интерлейкин-6, воспаление и иммунитет",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Высокая продукция IL-6 - провоспалительный профиль"),
                    "CG": ("moderate", "Умеренно повышенная продукция IL-6"),
                    "GC": ("moderate", "Умеренно повышенная продукция IL-6"),
                    "GG": ("normal", "Нормальная продукция IL-6"),
                }
            },
            "rs1800896": {
                "gene": "IL-10",
                "description": "Интерлейкин-10, противовоспалительный",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("low", "Низкая продукция IL-10 - сниженный противовоспалительный ответ"),
                    "AG": ("moderate", "Умеренная продукция IL-10"),
                    "GA": ("moderate", "Умеренная продукция IL-10"),
                    "GG": ("normal", "Высокая продукция IL-10 - хороший противовоспалительный ответ"),
                }
            },
            "rs16944": {
                "gene": "IL-1beta",
                "description": "Интерлейкин-1beta, воспаление",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Высокая продукция IL-1beta - провоспалительный профиль"),
                    "AG": ("moderate", "Умеренно повышенная продукция IL-1beta"),
                    "GA": ("moderate", "Умеренно повышенная продукция IL-1beta"),
                    "GG": ("normal", "Нормальная продукция IL-1beta"),
                }
            },
            "rs20541": {
                "gene": "IL-13",
                "description": "Интерлейкин-13, аллергия и астма",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Высокая продукция IL-13 - риск астмы и аллергии"),
                    "AG": ("moderate", "Повышенная продукция IL-13"),
                    "GA": ("moderate", "Повышенная продукция IL-13"),
                    "GG": ("normal", "Нормальная продукция IL-13"),
                }
            },
        }
    },

    "innate_immunity": {
        "name": "Врождённый иммунитет (TLR, комплемент)",
        "snps": {
            "rs5743708": {
                "gene": "TLR2",
                "description": "Toll-like рецептор 2, бактериальные инфекции",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Нарушение функции TLR2 - сниженная защита от бактерий"),
                    "AG": ("moderate", "Сниженная функция TLR2"),
                    "GA": ("moderate", "Сниженная функция TLR2"),
                    "GG": ("normal", "Нормальная функция TLR2"),
                }
            },
            "rs4986790": {
                "gene": "TLR4",
                "description": "Toll-like рецептор 4, грам-отрицательные бактерии",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("high", "Нарушение функции TLR4 - сниженный ответ на LPS"),
                    "AG": ("moderate", "Сниженная функция TLR4"),
                    "GA": ("moderate", "Сниженная функция TLR4"),
                    "AA": ("normal", "Нормальная функция TLR4"),
                }
            },
            "rs2230199": {
                "gene": "C3",
                "description": "Компонент комплемента C3",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("moderate", "C3F/F - повышенная активация комплемента, риск AMD"),
                    "CG": ("low", "C3S/F - умеренная активация"),
                    "GC": ("low", "C3S/F - умеренная активация"),
                    "CC": ("normal", "C3S/S - нормальная активация комплемента"),
                }
            },
        }
    },

    "infections": {
        "name": "Инфекционные заболевания",
        "snps": {
            "rs333": {
                "gene": "CCR5-delta32",
                "description": "Устойчивость к ВИЧ",
                "risk_allele": "D",
                "interpretation": {
                    "DD": ("protective", "CCR5-delta32 гомозигота - высокая устойчивость к ВИЧ-1"),
                    "DI": ("protective", "CCR5-delta32 гетерозигота - частичная защита от ВИЧ"),
                    "ID": ("protective", "CCR5-delta32 гетерозигота - частичная защита от ВИЧ"),
                    "--": ("protective", "Делеция - устойчивость к ВИЧ"),
                    "II": ("normal", "Нет делеции CCR5 - стандартная восприимчивость"),
                }
            },
            "rs12979860": {
                "gene": "IL28B (IFNL3)",
                "description": "Гепатит C - спонтанное излечение и ответ на терапию",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("protective", "Хороший ответ на терапию гепатита C, высокий шанс излечения"),
                    "CT": ("moderate", "Промежуточный ответ на терапию гепатита C"),
                    "TC": ("moderate", "Промежуточный ответ на терапию гепатита C"),
                    "TT": ("low", "Плохой ответ на терапию гепатита C"),
                }
            },
            "rs601338": {
                "gene": "FUT2",
                "description": "Норовирус, ротавирус - секреторный статус",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("protective", "Несекретор - устойчивость к норовирусу и ротавирусу"),
                    "AG": ("normal", "Секретор - стандартная восприимчивость"),
                    "GA": ("normal", "Секретор - стандартная восприимчивость"),
                    "GG": ("normal", "Секретор - стандартная восприимчивость к норовирусу"),
                }
            },
            "rs2814778": {
                "gene": "DARC (Duffy)",
                "description": "Малярия Plasmodium vivax",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("protective", "Duffy-отрицательный - устойчивость к P. vivax малярии"),
                    "CT": ("moderate", "Частичная защита от P. vivax"),
                    "TC": ("moderate", "Частичная защита от P. vivax"),
                    "TT": ("normal", "Duffy-положительный - восприимчивость к P. vivax"),
                }
            },
        }
    },

    "allergy": {
        "name": "Аллергия и атопия",
        "snps": {
            "rs7216389": {
                "gene": "ORMDL3",
                "description": "Детская астма",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Высокий риск детской астмы"),
                    "CT": ("moderate", "Повышенный риск астмы"),
                    "TC": ("moderate", "Повышенный риск астмы"),
                    "CC": ("normal", "Нормальный риск астмы"),
                }
            },
            "rs61816761": {
                "gene": "FLG (филаггрин)",
                "description": "Атопический дерматит, экзема",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Нарушение барьерной функции кожи - высокий риск экземы"),
                    "AG": ("moderate", "Носитель - повышенный риск атопического дерматита"),
                    "GA": ("moderate", "Носитель - повышенный риск атопического дерматита"),
                    "GG": ("normal", "Нормальная функция филаггрина"),
                }
            },
            "rs1801275": {
                "gene": "IL4RA",
                "description": "Рецептор IL-4, атопия и аллергия",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Повышенная чувствительность к IL-4 - риск атопии"),
                    "AG": ("moderate", "Умеренно повышенный риск аллергии"),
                    "GA": ("moderate", "Умеренно повышенный риск аллергии"),
                    "GG": ("normal", "Нормальная чувствительность к IL-4"),
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


def determine_celiac_risk(results):
    """Determine combined celiac disease risk from HLA-DQ2.5 and HLA-DQ8"""
    dq2 = None
    dq8 = None

    for r in results:
        if r['snp_id'] == 'rs2187668':
            dq2 = r
        elif r['snp_id'] == 'rs7454108':
            dq8 = r

    if not dq2 or not dq8:
        return None

    dq2_gt = dq2.get('genotype', '')
    dq8_gt = dq8.get('genotype', '')

    # Risk assessment
    dq2_risk = 'T' in dq2_gt if dq2_gt else False
    dq8_risk = 'C' in dq8_gt if dq8_gt else False

    if dq2_gt == 'TT':
        status = ('very_high', 'HLA-DQ2.5 гомозигота - очень высокий риск целиакии (>50%)')
    elif dq2_risk and dq8_risk:
        status = ('high', 'DQ2.5 + DQ8 - высокий риск целиакии')
    elif dq2_risk:
        status = ('moderate', 'HLA-DQ2.5 носитель - повышенный риск целиакии (~5-10%)')
    elif dq8_risk:
        status = ('low', 'HLA-DQ8 - небольшой риск целиакии (~2%)')
    else:
        status = ('normal', 'Низкий риск целиакии (<1%)')

    return {
        'dq2_genotype': dq2_gt,
        'dq8_genotype': dq8_gt,
        'status': status[0],
        'interpretation': status[1]
    }


def determine_inflammation_profile(results):
    """Determine overall inflammation profile from cytokine SNPs"""
    proinflammatory = 0
    antiinflammatory = 0
    total = 0

    cytokine_info = []

    for r in results:
        if not r['found']:
            continue
        total += 1

        if r['snp_id'] == 'rs1800629':  # TNF-alpha
            if 'A' in r['genotype']:
                proinflammatory += 1
                cytokine_info.append(f"TNF-alpha ({r['genotype']}): повышен")
            else:
                cytokine_info.append(f"TNF-alpha ({r['genotype']}): норма")

        elif r['snp_id'] == 'rs1800795':  # IL-6
            if 'C' in r['genotype']:
                proinflammatory += 1
                cytokine_info.append(f"IL-6 ({r['genotype']}): повышен")
            else:
                cytokine_info.append(f"IL-6 ({r['genotype']}): норма")

        elif r['snp_id'] == 'rs1800896':  # IL-10 (anti-inflammatory)
            if r['genotype'] == 'GG':
                antiinflammatory += 1
                cytokine_info.append(f"IL-10 ({r['genotype']}): высокий (защитный)")
            elif 'A' in r['genotype']:
                cytokine_info.append(f"IL-10 ({r['genotype']}): снижен")

        elif r['snp_id'] == 'rs16944':  # IL-1beta
            if 'A' in r['genotype']:
                proinflammatory += 1
                cytokine_info.append(f"IL-1beta ({r['genotype']}): повышен")
            else:
                cytokine_info.append(f"IL-1beta ({r['genotype']}): норма")

        elif r['snp_id'] == 'rs20541':  # IL-13
            if 'A' in r['genotype']:
                proinflammatory += 1
                cytokine_info.append(f"IL-13 ({r['genotype']}): повышен (аллергия)")
            else:
                cytokine_info.append(f"IL-13 ({r['genotype']}): норма")

    if total == 0:
        return None

    if proinflammatory >= 3:
        profile = ('high_inflammation', 'Провоспалительный профиль - рекомендуется противовоспалительная диета')
    elif proinflammatory >= 2:
        profile = ('moderate_inflammation', 'Умеренно провоспалительный профиль')
    elif antiinflammatory >= 1 and proinflammatory <= 1:
        profile = ('balanced', 'Сбалансированный воспалительный профиль')
    else:
        profile = ('normal', 'Нормальный воспалительный профиль')

    return {
        'proinflammatory_count': proinflammatory,
        'antiinflammatory_count': antiinflammatory,
        'total_analyzed': total,
        'profile': profile[0],
        'interpretation': profile[1],
        'details': cytokine_info
    }


def generate_category_report(category, results, genome):
    """Generate report for a category"""
    cat_info = IMMUNITY_SNPS[category]

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
        report.append("### Сводка по рискам\n")
        risk_emoji = {
            'high': '🔴',
            'very_high': '🔴🔴',
            'moderate': '🟡',
            'low': '🟢',
            'normal': '✅',
            'protective': '🛡️',
            'info': 'ℹ️'
        }
        for risk, count in sorted(risk_counts.items()):
            emoji = risk_emoji.get(risk, '•')
            report.append(f"- {emoji} {risk}: {count}")

    report.append("\n### Детальные результаты\n")
    report.append("| SNP | Ген | Генотип | Риск | Интерпретация |")
    report.append("|-----|-----|---------|------|---------------|")

    for r in results:
        if r['found']:
            risk_label = r['risk_level'] or 'н/д'
            interp = r['interpretation'] or 'Нет данных'
            report.append(f"| {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {risk_label} | {interp} |")
        else:
            report.append(f"| {r['snp_id']} | {r['gene']} | - | - | Не найден в геноме |")

    return '\n'.join(report)


def generate_summary_report(all_results, genome):
    """Generate overall immunity summary report"""
    report = []
    report.append("# Анализ иммунитета")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("\n---\n")

    report.append("## Важные предупреждения\n")
    report.append("1. **Это НЕ медицинский диагноз** - только информационный анализ")
    report.append("2. **Наличие риск-аллеля НЕ равно заболеванию** - пенетрантность варьируется")
    report.append("3. **Иммунитет зависит от многих факторов** - гены + образ жизни + среда")
    report.append("4. **Для медицинских решений** - консультация иммунолога/аллерголога обязательна\n")

    report.append("---\n")

    # Collect findings by risk level
    high_risk = []
    moderate_risk = []
    protective = []

    for category, results in all_results.items():
        for r in results:
            if r['risk_level'] in ['high', 'very_high']:
                high_risk.append((category, r))
            elif r['risk_level'] == 'moderate':
                moderate_risk.append((category, r))
            elif r['risk_level'] == 'protective':
                protective.append((category, r))

    if high_risk:
        report.append("## Маркеры повышенного риска\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in high_risk:
            cat_name = IMMUNITY_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if moderate_risk:
        report.append("## Маркеры умеренного риска\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in moderate_risk:
            cat_name = IMMUNITY_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if protective:
        report.append("## Защитные варианты\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in protective:
            cat_name = IMMUNITY_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    # Special analyses
    report.append("---\n")
    report.append("## Специальные анализы\n")

    # Celiac disease risk
    hla_results = all_results.get('hla_system', [])
    celiac = determine_celiac_risk(hla_results)
    if celiac:
        report.append("### Риск целиакии (HLA-DQ2.5/DQ8)\n")
        report.append(f"- HLA-DQ2.5 (rs2187668): {celiac['dq2_genotype']}")
        report.append(f"- HLA-DQ8 (rs7454108): {celiac['dq8_genotype']}")
        report.append(f"- **Статус: {celiac['status']}**")
        report.append(f"- {celiac['interpretation']}\n")

    # Inflammation profile
    cytokine_results = all_results.get('cytokines', [])
    inflammation = determine_inflammation_profile(cytokine_results)
    if inflammation:
        report.append("### Воспалительный профиль (цитокины)\n")
        report.append(f"- Провоспалительных маркеров: {inflammation['proinflammatory_count']}")
        report.append(f"- Противовоспалительных: {inflammation['antiinflammatory_count']}")
        report.append(f"- **Профиль: {inflammation['profile']}**")
        report.append(f"- {inflammation['interpretation']}")
        report.append("\nДетали:")
        for detail in inflammation['details']:
            report.append(f"  - {detail}")
        report.append("")

    # Infection resistance summary
    report.append("### Устойчивость к инфекциям\n")
    infection_results = all_results.get('infections', [])
    for r in infection_results:
        if r['found']:
            status = "защитный" if r['risk_level'] == 'protective' else r['risk_level'] or 'н/д'
            report.append(f"- **{r['gene']}** ({r['genotype']}): {r['interpretation']}")

    report.append("\n---\n")
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
    report.append("## Рекомендации\n")
    report.append("### При повышенных провоспалительных маркерах:")
    report.append("- Противовоспалительная диета (омега-3, куркума, имбирь)")
    report.append("- Контроль уровня витамина D")
    report.append("- Регулярная физическая активность")
    report.append("- Управление стрессом\n")
    report.append("### При риске аутоиммунных заболеваний:")
    report.append("- Регулярный мониторинг аутоантител")
    report.append("- Избегание триггеров (стресс, инфекции)")
    report.append("- Консультация ревматолога/иммунолога\n")
    report.append("### При риске аллергии/атопии:")
    report.append("- Укрепление кожного барьера (эмоленты)")
    report.append("- Идентификация и избегание аллергенов")
    report.append("- Пробиотики для иммунной модуляции")

    return '\n'.join(report)


def main():
    print("=" * 60)
    print("АНАЛИЗ ИММУНИТЕТА ПО ГЕНОМУ 23andMe")
    print("=" * 60)

    print("\n[1/4] Загрузка генома...")
    genome = load_genome()
    print(f"      Загружено {len(genome)} SNP")

    print("\n[2/4] Анализ маркеров по категориям...")
    all_results = {}

    for category, cat_info in IMMUNITY_SNPS.items():
        print(f"      -> {cat_info['name']}...")
        results = []
        for snp_id, snp_info in cat_info['snps'].items():
            result = analyze_snp(snp_id, snp_info, genome)
            results.append(result)
        all_results[category] = results

        # Count found
        found = sum(1 for r in results if r['found'])
        print(f"        Найдено: {found}/{len(results)}")

    print("\n[3/4] Генерация отчётов по категориям...")
    for category, results in all_results.items():
        report = generate_category_report(category, results, genome)
        report_dir = f"{REPORTS_PATH}/{category}"
        os.makedirs(report_dir, exist_ok=True)
        report_path = f"{report_dir}/report.md"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"      -> {report_path}")

    print("\n[4/4] Генерация сводного отчёта...")
    summary = generate_summary_report(all_results, genome)
    summary_dir = f"{REPORTS_PATH}/immunity"
    os.makedirs(summary_dir, exist_ok=True)
    summary_path = f"{summary_dir}/report.md"
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(summary)
    print(f"      -> {summary_path}")

    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)

    # Print key findings to console
    print("\nКЛЮЧЕВЫЕ НАХОДКИ:\n")

    for category, results in all_results.items():
        high_risk = [r for r in results if r['risk_level'] in ['high', 'very_high']]
        protective_found = [r for r in results if r['risk_level'] == 'protective']

        if high_risk:
            print(f"  {IMMUNITY_SNPS[category]['name']}:")
            for r in high_risk:
                print(f"    * {r['gene']} ({r['genotype']}): {r['interpretation']}")
            print()

        if protective_found:
            print(f"  {IMMUNITY_SNPS[category]['name']} (защитные):")
            for r in protective_found:
                print(f"    + {r['gene']} ({r['genotype']}): {r['interpretation']}")
            print()


if __name__ == "__main__":
    main()
