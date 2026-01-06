#!/usr/bin/env python3
"""
Reproductive SNP Analysis Script
Analyzes reproductive health genetic markers from 23andMe data
"""

import os
from collections import defaultdict
from datetime import datetime

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Organized by reproductive health category
# =============================================================================

REPRODUCTIVE_SNPS = {
    "female_ovarian": {
        "name": "Овариальный резерв (женщины)",
        "snps": {
            "rs16991615": {
                "gene": "MCM8",
                "description": "Время менопаузы",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Ранняя менопауза - значительно раньше среднего"),
                    "AG": ("moderate", "Менопауза немного раньше среднего (~1 год)"),
                    "GA": ("moderate", "Менопауза немного раньше среднего (~1 год)"),
                    "GG": ("normal", "Нормальное время менопаузы"),
                }
            },
            "rs244715": {
                "gene": "BRSK1",
                "description": "Овариальный резерв",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "Возможно сниженный овариальный резерв"),
                    "CT": ("low", "Незначительное влияние на резерв"),
                    "TC": ("low", "Незначительное влияние на резерв"),
                    "TT": ("normal", "Нормальный овариальный резерв"),
                }
            },
            "rs2303369": {
                "gene": "FNDC4",
                "description": "Функция яичников",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Возможно раннее снижение функции яичников"),
                    "CT": ("low", "Незначительное влияние"),
                    "TC": ("low", "Незначительное влияние"),
                    "CC": ("normal", "Нормальная функция яичников"),
                }
            },
        }
    },

    "pcos": {
        "name": "Синдром поликистозных яичников (СПКЯ)",
        "snps": {
            "rs13405728": {
                "gene": "LHCGR",
                "description": "Рецептор ЛГ/ХГ - риск СПКЯ",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Повышенный риск СПКЯ"),
                    "AG": ("moderate", "Умеренно повышенный риск СПКЯ"),
                    "GA": ("moderate", "Умеренно повышенный риск СПКЯ"),
                    "GG": ("normal", "Обычный риск"),
                }
            },
            "rs13429458": {
                "gene": "LHCGR",
                "description": "Рецептор ЛГ/ХГ - чувствительность",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Повышенный риск СПКЯ"),
                    "AG": ("moderate", "Умеренно повышенный риск"),
                    "GA": ("moderate", "Умеренно повышенный риск"),
                    "GG": ("normal", "Обычный риск"),
                }
            },
            "rs2479106": {
                "gene": "DENND1A",
                "description": "Ключевой ген СПКЯ",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Значительно повышенный риск СПКЯ"),
                    "AG": ("moderate", "Повышенный риск СПКЯ"),
                    "GA": ("moderate", "Повышенный риск СПКЯ"),
                    "GG": ("normal", "Обычный риск"),
                }
            },
            "rs10818854": {
                "gene": "DENND1A",
                "description": "Гиперандрогения при СПКЯ",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Повышенный риск гиперандрогении"),
                    "AG": ("moderate", "Умеренный риск гиперандрогении"),
                    "GA": ("moderate", "Умеренный риск гиперандрогении"),
                    "GG": ("normal", "Обычный уровень андрогенов"),
                }
            },
            "rs1801282": {
                "gene": "PPARG",
                "description": "Инсулинорезистентность при СПКЯ",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("info", "Pro/Pro - стандартный вариант"),
                    "CG": ("protective", "Pro/Ala - защита от инсулинорезистентности"),
                    "GC": ("protective", "Pro/Ala - защита от инсулинорезистентности"),
                    "GG": ("protective", "Ala/Ala - защитный эффект"),
                }
            },
        }
    },

    "pregnancy_risks": {
        "name": "Риски беременности",
        "snps": {
            "rs6025": {
                "gene": "F5 Leiden",
                "description": "Тромбофилия (риск тромбозов при беременности)",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Фактор V Лейден гомозигота - риск тромбозов 50x! Требуется антикоагуляция"),
                    "AG": ("high", "Фактор V Лейден гетерозигота - риск тромбозов 5x, наблюдение обязательно"),
                    "GA": ("high", "Фактор V Лейден гетерозигота - риск тромбозов 5x, наблюдение обязательно"),
                    "AC": ("high", "Носитель Factor V Leiden - повышенный риск тромбозов"),
                    "CA": ("high", "Носитель Factor V Leiden - повышенный риск тромбозов"),
                    "GG": ("normal", "Нет мутации Factor V Leiden"),
                    "CC": ("normal", "Нет мутации Factor V Leiden"),
                    "CG": ("normal", "Нет мутации Factor V Leiden"),
                    "GC": ("normal", "Нет мутации Factor V Leiden"),
                }
            },
            "rs1799963": {
                "gene": "F2 Prothrombin",
                "description": "Протромбин G20210A - риск венозных тромбозов",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Гомозигота - очень высокий риск тромбозов при беременности"),
                    "AG": ("high", "Носитель - повышенный риск тромбозов, требуется наблюдение"),
                    "GA": ("high", "Носитель - повышенный риск тромбозов, требуется наблюдение"),
                    "GG": ("normal", "Нет мутации протромбина"),
                }
            },
            "rs1801133": {
                "gene": "MTHFR C677T",
                "description": "Фолатный метаболизм - риск дефектов нервной трубки",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Сниженная активность MTHFR (~30%) - требуется метилфолат"),
                    "CT": ("moderate", "Умеренно сниженная активность (~65%) - рекомендован метилфолат"),
                    "TC": ("moderate", "Умеренно сниженная активность (~65%) - рекомендован метилфолат"),
                    "AG": ("moderate", "Гетерозигота - умеренно сниженная активность"),
                    "GA": ("moderate", "Гетерозигота - умеренно сниженная активность"),
                    "CC": ("normal", "Нормальная активность MTHFR"),
                    "AA": ("normal", "Нормальная активность MTHFR"),
                }
            },
            "rs1799889": {
                "gene": "SERPINE1 (PAI-1)",
                "description": "Риск преэклампсии и невынашивания",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("high", "4G/4G - повышенный риск преэклампсии и тромбозов"),
                    "AG": ("moderate", "4G/5G - умеренный риск"),
                    "GA": ("moderate", "4G/5G - умеренный риск"),
                    "AA": ("normal", "5G/5G - нормальный риск"),
                }
            },
        }
    },

    "gestational_diabetes": {
        "name": "Гестационный диабет",
        "snps": {
            "rs7903146": {
                "gene": "TCF7L2",
                "description": "Главный ген риска диабета 2 типа и гестационного диабета",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Значительно повышен риск ГСД (~80% выше)"),
                    "CT": ("moderate", "Умеренно повышен риск ГСД (~40% выше)"),
                    "TC": ("moderate", "Умеренно повышен риск ГСД (~40% выше)"),
                    "CC": ("normal", "Обычный риск ГСД"),
                }
            },
            "rs10830963": {
                "gene": "MTNR1B",
                "description": "Рецептор мелатонина - инсулиновый ответ",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("high", "Повышенный риск ГСД, нарушение секреции инсулина"),
                    "CG": ("moderate", "Умеренный риск ГСД"),
                    "GC": ("moderate", "Умеренный риск ГСД"),
                    "CC": ("normal", "Обычный риск"),
                }
            },
        }
    },

    "estrogen_metabolism": {
        "name": "Метаболизм эстрогенов",
        "snps": {
            "rs10046": {
                "gene": "CYP19A1",
                "description": "Ароматаза - синтез эстрогенов",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("info", "Повышенная активность ароматазы, выше уровень эстрогенов"),
                    "CT": ("info", "Средняя активность ароматазы"),
                    "TC": ("info", "Средняя активность ароматазы"),
                    "CC": ("info", "Сниженная активность ароматазы, ниже эстрогены"),
                }
            },
            "rs1056836": {
                "gene": "CYP1B1",
                "description": "4-OH эстрогены (потенциально генотоксичные)",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("moderate", "Val/Val - повышенное образование 4-OH эстрогенов, выше риск"),
                    "CG": ("low", "Leu/Val - умеренное образование"),
                    "GC": ("low", "Leu/Val - умеренное образование"),
                    "CC": ("normal", "Leu/Leu - нормальный метаболизм эстрогенов"),
                }
            },
            "rs4680": {
                "gene": "COMT",
                "description": "Метилирование эстрогенов - детоксикация",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("moderate", "Met/Met - медленное метилирование эстрогенов, дольше в организме"),
                    "AG": ("info", "Val/Met - среднее метилирование"),
                    "GA": ("info", "Val/Met - среднее метилирование"),
                    "GG": ("info", "Val/Val - быстрое метилирование эстрогенов"),
                }
            },
        }
    },

    "male_testosterone": {
        "name": "Тестостерон (мужчины)",
        "snps": {
            "rs6258": {
                "gene": "SHBG",
                "description": "Связывающий глобулин - уровень свободного тестостерона",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("info", "Снижен SHBG - выше свободный тестостерон"),
                    "CT": ("info", "Средний уровень SHBG"),
                    "TC": ("info", "Средний уровень SHBG"),
                    "CC": ("info", "Нормальный SHBG"),
                }
            },
            "rs12150660": {
                "gene": "SHBG",
                "description": "Уровень SHBG и тестостерона",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("info", "Низкий SHBG - выше биодоступный тестостерон"),
                    "GT": ("info", "Средний SHBG"),
                    "TG": ("info", "Средний SHBG"),
                    "GG": ("info", "Высокий SHBG - ниже свободный тестостерон"),
                }
            },
            "rs727428": {
                "gene": "SHBG",
                "description": "Уровень связывающего глобулина",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("info", "Низкий SHBG"),
                    "AG": ("info", "Средний SHBG"),
                    "GA": ("info", "Средний SHBG"),
                    "GG": ("info", "Высокий SHBG"),
                }
            },
            "rs10046": {
                "gene": "CYP19A1",
                "description": "Ароматаза - конверсия тестостерона в эстрогены",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Высокая ароматаза - больше конверсия в эстрогены"),
                    "CT": ("info", "Средняя активность ароматазы"),
                    "TC": ("info", "Средняя активность ароматазы"),
                    "CC": ("info", "Низкая ароматаза - меньше конверсия"),
                }
            },
        }
    },

    "fertility_male": {
        "name": "Мужская фертильность",
        "snps": {
            "rs5934505": {
                "gene": "FSHB",
                "description": "Фолликулостимулирующий гормон - сперматогенез",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Сниженная продукция ФСГ - возможно влияние на сперматогенез"),
                    "CT": ("low", "Незначительное снижение ФСГ"),
                    "TC": ("low", "Незначительное снижение ФСГ"),
                    "CC": ("normal", "Нормальная продукция ФСГ"),
                }
            },
            "rs1801133": {
                "gene": "MTHFR C677T",
                "description": "Качество спермы и ДНК фрагментация",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Сниженное качество спермы, повышенная ДНК фрагментация"),
                    "CT": ("low", "Незначительное влияние на качество"),
                    "TC": ("low", "Незначительное влияние на качество"),
                    "AG": ("low", "Незначительное влияние на качество спермы"),
                    "GA": ("low", "Незначительное влияние на качество спермы"),
                    "CC": ("normal", "Нормальное качество спермы"),
                    "AA": ("normal", "Нормальное качество спермы"),
                }
            },
            "rs4880": {
                "gene": "SOD2",
                "description": "Антиоксидант - защита спермы от окислительного стресса",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Ala/Ala - сниженная защита от оксидативного стресса"),
                    "CT": ("low", "Val/Ala - средняя защита"),
                    "TC": ("low", "Val/Ala - средняя защита"),
                    "CC": ("normal", "Val/Val - нормальная антиоксидантная защита"),
                }
            },
        }
    },

    "brca": {
        "name": "BRCA онкориски (репродуктивные органы)",
        "snps": {
            "rs80357906": {
                "gene": "BRCA1 (5382insC)",
                "description": "Мутация BRCA1 - рак молочной железы и яичников",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Носитель патогенной мутации BRCA1 - требуется генетическое консультирование!"),
                    "CT": ("high", "Носитель патогенной мутации BRCA1 - требуется генетическое консультирование!"),
                    "TC": ("high", "Носитель патогенной мутации BRCA1 - требуется генетическое консультирование!"),
                    "TT": ("normal", "Нет мутации 5382insC"),
                    "--": ("normal", "Нет мутации"),
                }
            },
            "rs28897672": {
                "gene": "BRCA1 (C61G)",
                "description": "Мутация BRCA1 - высокопенетрантная",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Носитель патогенной мутации BRCA1!"),
                    "GT": ("high", "Носитель патогенной мутации BRCA1!"),
                    "TG": ("high", "Носитель патогенной мутации BRCA1!"),
                    "GG": ("normal", "Нет мутации C61G"),
                    "--": ("normal", "Нет мутации"),
                }
            },
            "rs80359550": {
                "gene": "BRCA2 (6174delT)",
                "description": "Мутация BRCA2 - рак молочной железы, яичников, простаты",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Носитель патогенной мутации BRCA2!"),
                    "AT": ("high", "Носитель патогенной мутации BRCA2!"),
                    "TA": ("high", "Носитель патогенной мутации BRCA2!"),
                    "AA": ("normal", "Нет мутации 6174delT"),
                    "--": ("normal", "Нет мутации"),
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


def determine_thrombophilia_status(results):
    """Determine combined thrombophilia risk from Factor V and Prothrombin"""
    factor_v = None
    prothrombin = None

    for r in results:
        if r['snp_id'] == 'rs6025':
            factor_v = r
        elif r['snp_id'] == 'rs1799963':
            prothrombin = r

    if not factor_v or not prothrombin:
        return None

    fv_risk = factor_v.get('risk_level', 'normal')
    pt_risk = prothrombin.get('risk_level', 'normal')

    # Combined risk assessment
    if fv_risk == 'high' and pt_risk == 'high':
        status = ('very_high', 'Комбинированная тромбофилия - ОЧЕНЬ высокий риск! Обязательна антикоагуляция при беременности')
    elif fv_risk == 'high' or pt_risk == 'high':
        if factor_v.get('genotype', '') in ['AA']:
            status = ('very_high', 'Гомозигота Factor V Leiden - риск тромбозов 50x')
        else:
            status = ('high', 'Носительство тромбофилии - требуется наблюдение гематолога при беременности')
    else:
        status = ('normal', 'Нет наследственной тромбофилии')

    return {
        'factor_v': factor_v.get('genotype', 'N/A'),
        'prothrombin': prothrombin.get('genotype', 'N/A'),
        'status': status[0],
        'interpretation': status[1]
    }


def determine_estrogen_risk(results):
    """Determine estrogen metabolism risk profile"""
    cyp1b1 = None
    comt = None

    for r in results:
        if r['snp_id'] == 'rs1056836':
            cyp1b1 = r
        elif r['snp_id'] == 'rs4680':
            comt = r

    if not cyp1b1 or not comt:
        return None

    # CYP1B1 GG = high 4-OH, COMT AA = slow methylation = worst combination
    cyp1b1_high = cyp1b1.get('genotype', '') in ['GG', 'CG', 'GC']
    comt_slow = comt.get('genotype', '') in ['AA', 'AG', 'GA']

    if cyp1b1.get('genotype', '') == 'GG' and comt.get('genotype', '') == 'AA':
        status = ('high', 'Неблагоприятный профиль: высокие 4-OH эстрогены + медленное выведение')
    elif cyp1b1_high and comt_slow:
        status = ('moderate', 'Умеренный риск: повышенные 4-OH эстрогены, замедленное метилирование')
    elif cyp1b1_high or comt_slow:
        status = ('low', 'Один неблагоприятный фактор - следить за уровнем эстрогенов')
    else:
        status = ('normal', 'Благоприятный профиль метаболизма эстрогенов')

    return {
        'cyp1b1': cyp1b1.get('genotype', 'N/A'),
        'comt': comt.get('genotype', 'N/A'),
        'status': status[0],
        'interpretation': status[1]
    }


def generate_category_report(category, results, genome):
    """Generate report for a category"""
    cat_info = REPRODUCTIVE_SNPS[category]

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

    # Special sections based on category
    if category == 'pregnancy_risks':
        report.append("\n### Статус тромбофилии\n")
        thrombo = determine_thrombophilia_status(results)
        if thrombo:
            report.append(f"- Factor V Leiden (rs6025): {thrombo['factor_v']}")
            report.append(f"- Prothrombin G20210A (rs1799963): {thrombo['prothrombin']}")
            report.append(f"- **Статус: {thrombo['status']}**")
            report.append(f"- {thrombo['interpretation']}")

    if category == 'estrogen_metabolism':
        report.append("\n### Профиль метаболизма эстрогенов\n")
        estrogen = determine_estrogen_risk(results)
        if estrogen:
            report.append(f"- CYP1B1 (rs1056836): {estrogen['cyp1b1']} - образование 4-OH эстрогенов")
            report.append(f"- COMT (rs4680): {estrogen['comt']} - метилирование/выведение")
            report.append(f"- **Профиль: {estrogen['status']}**")
            report.append(f"- {estrogen['interpretation']}")

    return '\n'.join(report)


def generate_summary_report(all_results, genome):
    """Generate overall reproductive health summary report"""
    report = []
    report.append("# Сводный отчёт по репродуктивному здоровью")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("\n---\n")

    report.append("## Важные предупреждения\n")
    report.append("1. **Это НЕ медицинский диагноз** - только информационный анализ")
    report.append("2. **Генетика определяет предрасположенность**, а не судьбу")
    report.append("3. **Многие SNP имеют разную значимость** для мужчин и женщин")
    report.append("4. **Для планирования беременности** - консультация генетика обязательна")
    report.append("5. **BRCA мутации** требуют подтверждения клиническим тестированием\n")

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
        report.append("## 🔴 Маркеры повышенного риска\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in high_risk:
            cat_name = REPRODUCTIVE_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if moderate_risk:
        report.append("## 🟡 Маркеры умеренного риска\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in moderate_risk:
            cat_name = REPRODUCTIVE_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if protective:
        report.append("## 🛡️ Защитные варианты\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in protective:
            cat_name = REPRODUCTIVE_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    # Special analyses
    report.append("---\n")
    report.append("## Специальные анализы\n")

    # Thrombophilia status
    pregnancy_results = all_results.get('pregnancy_risks', [])
    thrombo = determine_thrombophilia_status(pregnancy_results)
    if thrombo:
        report.append("### Тромбофилия (риск при беременности)\n")
        report.append(f"- Factor V Leiden: {thrombo['factor_v']}")
        report.append(f"- Prothrombin: {thrombo['prothrombin']}")
        report.append(f"- **Статус: {thrombo['status']}**")
        report.append(f"- {thrombo['interpretation']}\n")

    # Estrogen metabolism
    estrogen_results = all_results.get('estrogen_metabolism', [])
    estrogen = determine_estrogen_risk(estrogen_results)
    if estrogen:
        report.append("### Метаболизм эстрогенов\n")
        report.append(f"- CYP1B1 (4-OH эстрогены): {estrogen['cyp1b1']}")
        report.append(f"- COMT (метилирование): {estrogen['comt']}")
        report.append(f"- **Профиль: {estrogen['status']}**")
        report.append(f"- {estrogen['interpretation']}\n")

    # BRCA summary
    report.append("### BRCA онкогены\n")
    brca_results = all_results.get('brca', [])
    brca_risk = False
    for r in brca_results:
        if r['found'] and r['risk_level'] == 'high':
            brca_risk = True
            report.append(f"- **ВНИМАНИЕ:** {r['gene']} ({r['genotype']}) - {r['interpretation']}")

    if not brca_risk:
        found_brca = [r for r in brca_results if r['found']]
        if found_brca:
            report.append("- Патогенные мутации BRCA1/BRCA2 не обнаружены")
        else:
            report.append("- Маркеры BRCA не найдены в геноме (требуется расширенное тестирование)")

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

    return '\n'.join(report)


def main():
    print("=" * 60)
    print("АНАЛИЗ РЕПРОДУКТИВНОГО ЗДОРОВЬЯ ПО ГЕНОМУ 23andMe")
    print("=" * 60)

    print("\n[1/4] Загрузка генома...")
    genome = load_genome()
    print(f"      Загружено {len(genome)} SNP")

    print("\n[2/4] Анализ маркеров по категориям...")
    all_results = {}

    for category, cat_info in REPRODUCTIVE_SNPS.items():
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
        category_dir = f"{REPORTS_PATH}/{category}"
        if not os.path.exists(category_dir):
            os.makedirs(category_dir)
        report_path = f"{category_dir}/report.md"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"      -> {report_path}")

    print("\n[4/4] Генерация сводного отчёта...")
    summary = generate_summary_report(all_results, genome)
    # Ensure reproductive directory exists
    repro_dir = f"{REPORTS_PATH}/reproductive"
    if not os.path.exists(repro_dir):
        os.makedirs(repro_dir)
    summary_path = f"{repro_dir}/report.md"
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
        if high_risk:
            print(f"  {REPRODUCTIVE_SNPS[category]['name']}:")
            for r in high_risk:
                print(f"    * {r['gene']} ({r['genotype']}): {r['interpretation']}")
            print()


if __name__ == "__main__":
    main()
