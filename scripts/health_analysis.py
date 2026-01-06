#!/usr/bin/env python3
"""
Health SNP Analysis Script
Analyzes health-related genetic markers from 23andMe data
"""

import os
from collections import defaultdict
from datetime import datetime

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Organized by health category
# =============================================================================

HEALTH_SNPS = {
    "cardiovascular": {
        "name": "Сердечно-сосудистые заболевания",
        "snps": {
            "rs10757274": {
                "gene": "9p21",
                "description": "ИБС, инфаркт миокарда",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("high", "Повышенный риск ИБС (гомозигота)"),
                    "AG": ("moderate", "Умеренно повышенный риск ИБС (гетерозигота)"),
                    "AA": ("normal", "Нормальный риск"),
                }
            },
            "rs1333049": {
                "gene": "9p21",
                "description": "ИБС, коронарная болезнь сердца",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Повышенный риск ИБС"),
                    "CG": ("moderate", "Умеренно повышенный риск"),
                    "GG": ("normal", "Нормальный риск"),
                }
            },
            "rs1801133": {
                "gene": "MTHFR C677T",
                "description": "Гомоцистеин, риск тромбов",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Сниженная активность MTHFR (~30%), повышен гомоцистеин"),
                    "CT": ("moderate", "Умеренно сниженная активность (~65%)"),
                    "AG": ("moderate", "Гетерозигота C677T - умеренно сниженная активность (~65%)"),
                    "CC": ("normal", "Нормальная активность MTHFR"),
                    "AA": ("normal", "Нормальная активность MTHFR"),
                }
            },
            "rs1801131": {
                "gene": "MTHFR A1298C",
                "description": "Метаболизм фолатов",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "Сниженная активность MTHFR"),
                    "AC": ("low", "Незначительное снижение"),
                    "GT": ("normal", "Нормальная активность A1298C"),
                    "TT": ("normal", "Нормальная активность"),
                    "AA": ("normal", "Нормальная активность"),
                }
            },
            "rs6025": {
                "gene": "Factor V Leiden",
                "description": "Тромбофилия",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Фактор V Лейден - высокий риск тромбозов!"),
                    "AG": ("high", "Носитель Factor V Leiden - повышенный риск тромбозов"),
                    "AC": ("high", "Носитель Factor V Leiden - повышенный риск тромбозов"),
                    "GG": ("normal", "Нет мутации Factor V Leiden"),
                    "CC": ("normal", "Нет мутации Factor V Leiden"),
                    "CG": ("normal", "Нет мутации Factor V Leiden"),
                }
            },
            "rs1799963": {
                "gene": "Prothrombin G20210A",
                "description": "Тромбозы венозные",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Гомозигота - высокий риск тромбозов"),
                    "AG": ("high", "Носитель - повышенный риск тромбозов"),
                    "GG": ("normal", "Нет мутации протромбина"),
                }
            },
            "rs1800562": {
                "gene": "HFE C282Y",
                "description": "Гемохроматоз (накопление железа)",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Гомозигота C282Y - высокий риск гемохроматоза"),
                    "AG": ("moderate", "Носитель C282Y"),
                    "GG": ("normal", "Нет мутации C282Y"),
                }
            },
        }
    },

    "oncology": {
        "name": "Онкологические маркеры",
        "snps": {
            "rs1042522": {
                "gene": "TP53 Arg72Pro",
                "description": "Общий риск рака",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "Pro/Pro - немного повышен риск некоторых раков"),
                    "CG": ("normal", "Arg/Pro - гетерозигота"),
                    "GG": ("normal", "Arg/Arg - стандартный вариант"),
                }
            },
            "rs1800566": {
                "gene": "NQO1",
                "description": "Детоксикация канцерогенов",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Нет активности NQO1 - сниженная детоксикация"),
                    "CT": ("moderate", "Сниженная активность NQO1"),
                    "CC": ("normal", "Нормальная детоксикация"),
                }
            },
            "rs1695": {
                "gene": "GSTP1",
                "description": "Детоксикация",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("moderate", "Сниженная детоксикация через GSTP1"),
                    "AG": ("low", "Немного сниженная активность"),
                    "AA": ("normal", "Нормальная активность GSTP1"),
                }
            },
            "rs1048943": {
                "gene": "CYP1A1",
                "description": "Метаболизм канцерогенов",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("moderate", "Повышенная активация канцерогенов"),
                    "AG": ("low", "Немного повышена"),
                    "AA": ("normal", "Нормальный метаболизм"),
                }
            },
        }
    },

    "neurology": {
        "name": "Неврология и психическое здоровье",
        "snps": {
            "rs429358": {
                "gene": "APOE (ε4 маркер 1)",
                "description": "Болезнь Альцгеймера",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Вероятно ε4/ε4 - значительно повышен риск"),
                    "CT": ("moderate", "Вероятно носитель ε4"),
                    "TT": ("normal", "Нет аллеля ε4"),
                }
            },
            "rs7412": {
                "gene": "APOE (ε маркер 2)",
                "description": "Болезнь Альцгеймера",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("info", "Используется для определения APOE генотипа"),
                    "CT": ("info", "Используется для определения APOE генотипа"),
                    "TT": ("info", "Используется для определения APOE генотипа"),
                }
            },
            "rs6265": {
                "gene": "BDNF Val66Met",
                "description": "Память, нейропластичность, депрессия",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("moderate", "Met/Met - сниженная секреция BDNF, риск депрессии"),
                    "TT": ("moderate", "Met/Met - сниженная секреция BDNF, риск депрессии"),
                    "AG": ("low", "Val/Met - немного сниженная секреция"),
                    "CT": ("low", "Val/Met - немного сниженная секреция"),
                    "GG": ("normal", "Val/Val - нормальная секреция BDNF"),
                    "CC": ("normal", "Val/Val - нормальная секреция BDNF"),
                }
            },
            "rs4680": {
                "gene": "COMT Val158Met",
                "description": "Стресс, тревожность, когнитивные функции",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("info", "Met/Met - 'Worrier' - лучше когнитивно, но выше тревожность"),
                    "AG": ("info", "Val/Met - сбалансированный тип"),
                    "GG": ("info", "Val/Val - 'Warrior' - устойчив к стрессу, хуже память"),
                }
            },
            "rs1800497": {
                "gene": "DRD2/ANKK1 Taq1A",
                "description": "Дофамин, зависимости",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("moderate", "A1/A1 - меньше D2 рецепторов, риск зависимостей"),
                    "AG": ("low", "A1/A2 - умеренно снижены D2 рецепторы"),
                    "GG": ("normal", "A2/A2 - нормальное количество D2 рецепторов"),
                }
            },
            "rs53576": {
                "gene": "OXTR",
                "description": "Окситоцин, эмпатия, социальность",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("info", "Ниже эмпатия и социальная чувствительность"),
                    "AG": ("info", "Средний уровень"),
                    "GG": ("info", "Выше эмпатия и социальная чувствительность"),
                }
            },
        }
    },

    "diabetes": {
        "name": "Диабет и метаболизм",
        "snps": {
            "rs7903146": {
                "gene": "TCF7L2",
                "description": "Диабет 2 типа (главный маркер)",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Значительно повышен риск СД2 (~80%)"),
                    "CT": ("moderate", "Повышен риск СД2 (~40%)"),
                    "CC": ("normal", "Нормальный риск"),
                }
            },
            "rs1801282": {
                "gene": "PPARG Pro12Ala",
                "description": "Инсулин, ожирение",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("normal", "Pro/Pro - стандартный вариант"),
                    "CG": ("protective", "Pro/Ala - защитный эффект против СД2"),
                    "GG": ("protective", "Ala/Ala - защитный эффект"),
                }
            },
            "rs5219": {
                "gene": "KCNJ11",
                "description": "Диабет 2 типа",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Повышен риск СД2"),
                    "CT": ("low", "Немного повышен риск"),
                    "CC": ("normal", "Нормальный риск"),
                }
            },
            "rs9939609": {
                "gene": "FTO",
                "description": "Ожирение (главный маркер)",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Повышен риск ожирения (+3кг в среднем)"),
                    "AT": ("moderate", "Умеренно повышен риск (+1.5кг)"),
                    "TT": ("normal", "Нормальный риск"),
                }
            },
            "rs17782313": {
                "gene": "MC4R",
                "description": "Аппетит, ожирение",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "Повышен аппетит и риск ожирения"),
                    "CT": ("low", "Немного повышен риск"),
                    "TT": ("normal", "Нормальный контроль аппетита"),
                }
            },
        }
    },

    "autoimmune": {
        "name": "Аутоиммунные заболевания",
        "snps": {
            "rs2187668": {
                "gene": "HLA-DQ2.5",
                "description": "Целиакия",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "HLA-DQ2.5 - высокий риск целиакии"),
                    "CT": ("moderate", "Носитель HLA-DQ2.5"),
                    "CC": ("normal", "Низкий риск целиакии"),
                }
            },
            "rs7454108": {
                "gene": "HLA-DQ8",
                "description": "Целиакия, диабет 1 типа",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "HLA-DQ8 - риск целиакии и СД1"),
                    "CT": ("low", "Носитель"),
                    "TT": ("normal", "Низкий риск"),
                }
            },
            "rs2476601": {
                "gene": "PTPN22",
                "description": "Аутоиммунные (RA, SLE, СД1)",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Высокий риск аутоиммунных заболеваний"),
                    "AG": ("moderate", "Повышен риск аутоиммунных"),
                    "GG": ("normal", "Нормальный риск"),
                }
            },
            "rs3135388": {
                "gene": "HLA-DRB1",
                "description": "Рассеянный склероз",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("moderate", "Повышен риск рассеянного склероза"),
                    "AG": ("low", "Немного повышен риск"),
                    "GG": ("normal", "Нормальный риск"),
                }
            },
        }
    },

    "pharmacogenomics": {
        "name": "Фармакогеномика",
        "snps": {
            "rs4244285": {
                "gene": "CYP2C19*2",
                "description": "Метаболизм: Клопидогрел, омепразол, антидепрессанты",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Плохой метаболизатор - клопидогрел НЕ эффективен!"),
                    "AG": ("moderate", "Промежуточный метаболизатор"),
                    "GG": ("normal", "Нормальный метаболизатор"),
                }
            },
            "rs4986893": {
                "gene": "CYP2C19*3",
                "description": "Метаболизм антидепрессантов",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Плохой метаболизатор"),
                    "AG": ("moderate", "Промежуточный"),
                    "GG": ("normal", "Нормальный метаболизатор"),
                }
            },
            "rs1065852": {
                "gene": "CYP2D6",
                "description": "Метаболизм: Кодеин, тамоксифен, антидепрессанты",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Плохой метаболизатор - кодеин не работает"),
                    "AG": ("moderate", "Промежуточный метаболизатор"),
                    "GG": ("normal", "Нормальный метаболизатор"),
                }
            },
            "rs9923231": {
                "gene": "VKORC1",
                "description": "Чувствительность к Варфарину",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Высокая чувствительность - нужна низкая доза варфарина"),
                    "CT": ("moderate", "Средняя чувствительность"),
                    "CC": ("normal", "Стандартная доза варфарина"),
                }
            },
            "rs1799853": {
                "gene": "CYP2C9*2",
                "description": "Метаболизм: Варфарин, НПВС",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("high", "Плохой метаболизатор - снизить дозу варфарина"),
                    "CT": ("moderate", "Промежуточный метаболизатор"),
                    "CC": ("normal", "Нормальный метаболизатор"),
                }
            },
            "rs1057910": {
                "gene": "CYP2C9*3",
                "description": "Метаболизм варфарина",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Плохой метаболизатор"),
                    "AC": ("moderate", "Промежуточный"),
                    "AA": ("normal", "Нормальный метаболизатор"),
                }
            },
            "rs4149056": {
                "gene": "SLCO1B1",
                "description": "Транспорт статинов (риск миопатии)",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("high", "Высокий риск миопатии от статинов!"),
                    "CT": ("moderate", "Повышен риск миопатии"),
                    "TT": ("normal", "Нормальный риск"),
                }
            },
            "rs1800460": {
                "gene": "TPMT",
                "description": "Метаболизм азатиоприна",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("high", "Плохой метаболизатор - токсичность азатиоприна!"),
                    "TT": ("high", "Плохой метаболизатор - токсичность азатиоприна!"),
                    "AG": ("moderate", "Промежуточный - снизить дозу"),
                    "CT": ("moderate", "Промежуточный - снизить дозу"),
                    "GG": ("normal", "Нормальный метаболизатор"),
                    "CC": ("normal", "Нормальный метаболизатор"),
                }
            },
        }
    },

    "other": {
        "name": "Другие важные маркеры",
        "snps": {
            "rs1800629": {
                "gene": "TNF-α",
                "description": "Воспаление",
                "risk_allele": "A",
                "interpretation": {
                    "AA": ("moderate", "Повышенная продукция TNF-α, хроническое воспаление"),
                    "AG": ("low", "Немного повышена"),
                    "GG": ("normal", "Нормальный уровень"),
                }
            },
            "rs1143634": {
                "gene": "IL-1β",
                "description": "Воспаление",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Повышенное воспаление"),
                    "CT": ("low", "Немного повышено"),
                    "AG": ("low", "Немного повышена продукция IL-1β"),
                    "CC": ("normal", "Нормальный уровень"),
                    "AA": ("normal", "Нормальный уровень"),
                }
            },
            "rs1800795": {
                "gene": "IL-6",
                "description": "Воспаление, иммунитет",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "Повышенная продукция IL-6"),
                    "CG": ("low", "Умеренно повышена"),
                    "GG": ("normal", "Нормальный уровень"),
                }
            },
            "rs7041": {
                "gene": "GC (VDBP)",
                "description": "Метаболизм витамина D",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Gc1F/1F - сниженный уровень витамина D"),
                    "GT": ("low", "Немного снижен"),
                    "AC": ("low", "Gc1S/1F - немного снижен витамин D"),
                    "GG": ("normal", "Gc1S/1S - нормальный уровень"),
                    "AA": ("normal", "Нормальный уровень"),
                    "CC": ("normal", "Нормальный уровень"),
                }
            },
            "rs2282679": {
                "gene": "GC",
                "description": "Уровень витамина D",
                "risk_allele": "C",
                "interpretation": {
                    "CC": ("moderate", "Сниженный витамин D"),
                    "AC": ("low", "Немного снижен"),
                    "GT": ("low", "Немного снижен уровень витамина D"),
                    "GG": ("normal", "Нормальный уровень"),
                    "TT": ("normal", "Нормальный уровень"),
                    "AA": ("normal", "Нормальный уровень"),
                }
            },
            "rs855791": {
                "gene": "TMPRSS6",
                "description": "Уровень железа",
                "risk_allele": "T",
                "interpretation": {
                    "TT": ("moderate", "Сниженный уровень железа"),
                    "CT": ("low", "Немного снижен"),
                    "AG": ("low", "Немного снижен уровень железа"),
                    "GG": ("normal", "Нормальный уровень железа"),
                    "AA": ("normal", "Нормальный уровень"),
                    "CC": ("normal", "Нормальный уровень"),
                }
            },
            "rs1799945": {
                "gene": "HFE H63D",
                "description": "Гемохроматоз (лёгкая форма)",
                "risk_allele": "G",
                "interpretation": {
                    "GG": ("moderate", "H63D гомозигота - лёгкий риск накопления железа"),
                    "CG": ("low", "Носитель H63D"),
                    "CC": ("normal", "Нет мутации H63D"),
                }
            },
        }
    }
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


def determine_apoe_genotype(genome):
    """Determine APOE genotype from rs429358 and rs7412"""
    rs429358 = genome.get('rs429358', {}).get('genotype', '')
    rs7412 = genome.get('rs7412', {}).get('genotype', '')

    # APOE determination table
    # rs429358 (C=ε4), rs7412 (T=ε2)
    apoe_table = {
        ('TT', 'CC'): ('ε2/ε2', 'protective', 'Защитный генотип - пониженный риск Альцгеймера'),
        ('TT', 'CT'): ('ε2/ε3', 'protective', 'Немного пониженный риск'),
        ('CT', 'CC'): ('ε2/ε4', 'moderate', 'Смешанный - один защитный, один рисковый аллель'),
        ('TT', 'TT'): ('ε3/ε3', 'normal', 'Наиболее распространённый генотип - обычный риск'),
        ('CT', 'CT'): ('ε3/ε4', 'high', 'Повышенный риск Альцгеймера (~3x)'),
        ('CC', 'TT'): ('ε4/ε4', 'very_high', 'Значительно повышенный риск Альцгеймера (~12x)'),
        ('CT', 'TT'): ('ε3/ε4', 'high', 'Повышенный риск Альцгеймера (~3x)'),
        ('CC', 'CT'): ('ε4/ε4 или ε3/ε4', 'high', 'Повышенный риск'),
    }

    # Normalize genotypes
    n429 = normalize_genotype(rs429358)
    n7412 = normalize_genotype(rs7412)

    for (g1, g2), (apoe, risk, desc) in apoe_table.items():
        if normalize_genotype(g1) == n429 and normalize_genotype(g2) == n7412:
            return {
                'rs429358': rs429358,
                'rs7412': rs7412,
                'apoe_genotype': apoe,
                'risk_level': risk,
                'interpretation': desc
            }

    return {
        'rs429358': rs429358,
        'rs7412': rs7412,
        'apoe_genotype': 'Не определён',
        'risk_level': 'unknown',
        'interpretation': f'Комбинация {rs429358}/{rs7412} не в таблице'
    }


def determine_mthfr_status(results):
    """Determine combined MTHFR status"""
    c677t = None
    a1298c = None

    for r in results:
        if r['snp_id'] == 'rs1801133':
            c677t = r['genotype']
        elif r['snp_id'] == 'rs1801131':
            a1298c = r['genotype']

    if not c677t or not a1298c:
        return None

    # Normalize genotypes - C677T: C>T mutation, A1298C: A>C mutation
    # Count risk alleles
    c677t_risk = c677t.count('T') if c677t else 0  # T is risk
    # For A1298C, the genotype might be reported as G/T (complement strand)
    # A>C on forward = T>G on reverse, so G is risk on reverse
    a1298c_risk = a1298c.count('C') + a1298c.count('G') if a1298c else 0

    # Check for homozygous risk
    is_c677t_homo = c677t in ['TT', 'AA']  # TT on forward, AA could be different strand
    is_c677t_hetero = 'T' in c677t and 'C' in c677t or 'A' in c677t and 'G' in c677t
    is_a1298c_homo = a1298c in ['CC', 'GG']
    is_a1298c_hetero = len(set(a1298c)) == 2 if a1298c else False

    # More lenient detection: AG for C677T = heterozygote (strand complement)
    if c677t == 'AG':
        is_c677t_hetero = True
    if a1298c == 'GT':
        is_a1298c_hetero = False  # GT likely means normal/normal on this assay

    if is_c677t_homo and is_a1298c_homo:
        status = ('severe', 'Значительное снижение активности MTHFR (~10-20%)')
    elif is_c677t_homo and is_a1298c_hetero:
        status = ('severe', 'Значительное снижение активности MTHFR')
    elif is_c677t_homo:
        status = ('moderate', 'C677T гомозигота - сниженная активность MTHFR (~30%)')
    elif is_c677t_hetero and is_a1298c_homo:
        status = ('moderate', 'Компаунд - умеренное снижение активности')
    elif is_c677t_hetero and is_a1298c_hetero:
        status = ('moderate', 'Компаунд гетерозигота - умеренное снижение')
    elif is_c677t_hetero:
        status = ('mild', 'C677T гетерозигота - незначительное снижение (~65%)')
    elif is_a1298c_homo:
        status = ('mild', 'A1298C гомозигота - незначительное снижение')
    else:
        status = ('normal', 'Нормальная активность MTHFR')

    return {
        'c677t': c677t,
        'a1298c': a1298c,
        'status': status[0],
        'interpretation': status[1]
    }


def generate_category_report(category, results, genome):
    """Generate report for a category"""
    cat_info = HEALTH_SNPS[category]

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

    # Special sections
    if category == 'neurology':
        report.append("\n### APOE генотип (риск Альцгеймера)\n")
        apoe = determine_apoe_genotype(genome)
        report.append(f"- rs429358: {apoe['rs429358']}")
        report.append(f"- rs7412: {apoe['rs7412']}")
        report.append(f"- **APOE генотип: {apoe['apoe_genotype']}**")
        report.append(f"- Риск: {apoe['risk_level']}")
        report.append(f"- {apoe['interpretation']}")

    if category == 'cardiovascular':
        mthfr = determine_mthfr_status(results)
        if mthfr:
            report.append("\n### MTHFR статус\n")
            report.append(f"- C677T (rs1801133): {mthfr['c677t']}")
            report.append(f"- A1298C (rs1801131): {mthfr['a1298c']}")
            report.append(f"- **Статус: {mthfr['status']}**")
            report.append(f"- {mthfr['interpretation']}")

    return '\n'.join(report)


def generate_summary_report(all_results, genome):
    """Generate overall summary report"""
    report = []
    report.append("# 📊 Сводный отчёт по здоровью")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("\n---\n")

    report.append("## ⚠️ Важные предупреждения\n")
    report.append("1. **Это НЕ медицинский диагноз** — только информационный анализ")
    report.append("2. **Наличие риск-аллеля ≠ заболевание** — пенетрантность варьируется")
    report.append("3. **Большинство болезней полигенные** — зависят от многих генов + среда")
    report.append("4. **Для медицинских решений** — консультация генетика обязательна\n")

    report.append("---\n")

    # Collect high-risk findings
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
            cat_name = HEALTH_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if moderate_risk:
        report.append("## 🟡 Маркеры умеренного риска\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in moderate_risk:
            cat_name = HEALTH_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    if protective:
        report.append("## 🛡️ Защитные варианты\n")
        report.append("| Категория | SNP | Ген | Генотип | Описание |")
        report.append("|-----------|-----|-----|---------|----------|")
        for cat, r in protective:
            cat_name = HEALTH_SNPS[cat]['name']
            report.append(f"| {cat_name} | {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {r['interpretation']} |")
        report.append("")

    # Special analyses
    report.append("---\n")
    report.append("## 🧬 Специальные анализы\n")

    # APOE
    apoe = determine_apoe_genotype(genome)
    report.append("### APOE (Болезнь Альцгеймера)\n")
    report.append(f"- **Генотип: {apoe['apoe_genotype']}**")
    report.append(f"- {apoe['interpretation']}\n")

    # MTHFR
    cardio_results = all_results.get('cardiovascular', [])
    mthfr = determine_mthfr_status(cardio_results)
    if mthfr:
        report.append("### MTHFR (Метаболизм фолатов)\n")
        report.append(f"- C677T: {mthfr['c677t']}, A1298C: {mthfr['a1298c']}")
        report.append(f"- **Статус: {mthfr['status']}**")
        report.append(f"- {mthfr['interpretation']}\n")

    # Pharmacogenomics summary
    report.append("### 💊 Фармакогеномика - ключевые находки\n")
    pharma = all_results.get('pharmacogenomics', [])
    important_drugs = []
    for r in pharma:
        if r['found'] and r['risk_level'] in ['high', 'moderate']:
            important_drugs.append(f"- **{r['gene']}** ({r['genotype']}): {r['interpretation']}")

    if important_drugs:
        report.extend(important_drugs)
    else:
        report.append("- Нет критических фармакогенетических находок")

    report.append("\n---\n")
    report.append("## 📈 Статистика анализа\n")

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
    print("АНАЛИЗ ЗДОРОВЬЯ ПО ГЕНОМУ 23andMe")
    print("=" * 60)

    print("\n[1/4] Загрузка генома...")
    genome = load_genome()
    print(f"      Загружено {len(genome)} SNP")

    print("\n[2/4] Анализ маркеров по категориям...")
    all_results = {}

    for category, cat_info in HEALTH_SNPS.items():
        print(f"      → {cat_info['name']}...")
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
        report_path = f"{REPORTS_PATH}/{category}/report.md"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"      → {report_path}")

    print("\n[4/4] Генерация сводного отчёта...")
    summary = generate_summary_report(all_results, genome)
    summary_path = f"{REPORTS_PATH}/health_summary.md"
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(summary)
    print(f"      → {summary_path}")

    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)

    # Print key findings to console
    print("\n🔑 КЛЮЧЕВЫЕ НАХОДКИ:\n")

    for category, results in all_results.items():
        high_risk = [r for r in results if r['risk_level'] in ['high', 'very_high']]
        if high_risk:
            print(f"⚠️  {HEALTH_SNPS[category]['name']}:")
            for r in high_risk:
                print(f"    • {r['gene']} ({r['genotype']}): {r['interpretation']}")
            print()


if __name__ == "__main__":
    main()
