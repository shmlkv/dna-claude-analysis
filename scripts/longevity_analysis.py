#!/usr/bin/env python3
"""
Longevity and Aging Analysis Script
Analyzes genetic markers associated with lifespan, aging, and anti-aging pathways from 23andMe data
"""

import os
from datetime import datetime
from collections import defaultdict

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Longevity and Aging Markers
# =============================================================================

LONGEVITY_SNPS = {
    "longevity_genes": {
        "name": "Longevity Genes",
        "name_ru": "Гены долголетия",
        "snps": {
            "rs2802292": {
                "gene": "FOXO3",
                "description": "Forkhead box O3 - Master longevity regulator",
                "description_ru": "Ключевой регулятор долголетия, контролирует стресс-резистентность",
                "protective_allele": "T",
                "interpretation": {
                    "TT": ("beneficial", "Longevity variant (+2.7 years)", "Вариант долголетия (+2.7 года к ожидаемой продолжительности жизни)", 2.7),
                    "GT": ("moderate", "Partial longevity benefit (+1.5 years)", "Частичный эффект долголетия (+1.5 года)", 1.5),
                    "GG": ("baseline", "Baseline lifespan", "Базовая продолжительность жизни", 0),
                }
            },
            "rs5882": {
                "gene": "CETP",
                "description": "Cholesteryl ester transfer protein - HDL metabolism",
                "description_ru": "Белок переноса холестерина, связан с долголетием ашкенази",
                "protective_allele": "G",
                "interpretation": {
                    "GG": ("beneficial", "Longevity variant - higher HDL", "Вариант долголетия - повышенный HDL холестерин", 1.5),
                    "AG": ("moderate", "Partial HDL benefit", "Умеренно повышенный HDL", 0.7),
                    "AA": ("baseline", "Standard HDL metabolism", "Стандартный метаболизм HDL", 0),
                }
            },
            "rs9536314": {
                "gene": "KLOTHO (KL)",
                "description": "Klotho - Anti-aging hormone, longevity factor",
                "description_ru": "Антивозрастной гормон, защищает от возрастных болезней",
                "protective_allele": "T",
                "interpretation": {
                    "TT": ("beneficial", "KL-VS variant - cognitive and longevity benefit", "Вариант KL-VS - защита когнитивных функций и долголетие", 1.8),
                    "GT": ("beneficial", "Heterozygote - optimal Klotho levels", "Гетерозигота - оптимальные уровни клото (лучший вариант)", 2.0),
                    "GG": ("baseline", "Standard Klotho", "Стандартные уровни клото", 0),
                }
            },
        }
    },

    "apoe": {
        "name": "APOE Genotype",
        "name_ru": "Генотип APOE",
        "snps": {
            "rs429358": {
                "gene": "APOE",
                "description": "APOE C112R - determines ε4 status",
                "description_ru": "APOE C112R - определяет статус ε4",
                "note": "C=ε4 allele component, T=ε2/ε3 component",
                "interpretation": {
                    "CC": ("risk", "Two ε4 components", "Два компонента ε4", -3.0),
                    "CT": ("moderate", "One ε4 component", "Один компонент ε4", -1.5),
                    "TT": ("baseline", "No ε4 component", "Нет компонента ε4", 0),
                }
            },
            "rs7412": {
                "gene": "APOE",
                "description": "APOE R158C - determines ε2 status",
                "description_ru": "APOE R158C - определяет статус ε2",
                "note": "T=ε2 allele (protective), C=ε3/ε4 component",
                "interpretation": {
                    "TT": ("beneficial", "Two ε2 components", "Два компонента ε2 (защитный)", 2.5),
                    "CT": ("beneficial", "One ε2 component", "Один компонент ε2 (защитный)", 1.2),
                    "CC": ("baseline", "No ε2 component", "Нет компонента ε2", 0),
                }
            },
        }
    },

    "telomeres": {
        "name": "Telomere Length",
        "name_ru": "Длина теломер",
        "snps": {
            "rs10936599": {
                "gene": "TERC",
                "description": "Telomerase RNA component - telomere maintenance",
                "description_ru": "РНК-компонент теломеразы - поддержание теломер",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Longer telomeres", "Более длинные теломеры", 1.0),
                    "CT": ("moderate", "Average telomere length", "Средняя длина теломер", 0.5),
                    "TT": ("risk", "Shorter telomeres", "Укороченные теломеры", -0.5),
                }
            },
            "rs2736100": {
                "gene": "TERT",
                "description": "Telomerase reverse transcriptase - telomere extension",
                "description_ru": "Обратная транскриптаза теломеразы - удлинение теломер",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Enhanced telomerase activity", "Повышенная активность теломеразы", 1.0),
                    "AC": ("moderate", "Average telomerase", "Средняя активность теломеразы", 0.5),
                    "AA": ("baseline", "Standard telomerase", "Стандартная активность", 0),
                }
            },
        }
    },

    "oxidative_stress": {
        "name": "Oxidative Stress Defense",
        "name_ru": "Защита от оксидативного стресса",
        "snps": {
            "rs4880": {
                "gene": "SOD2 (MnSOD)",
                "description": "Superoxide dismutase 2 - mitochondrial antioxidant",
                "description_ru": "Супероксиддисмутаза 2 - митохондриальный антиоксидант",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Better antioxidant protection (Ala/Ala)", "Лучшая антиоксидантная защита (Ala/Ala)", 1.0),
                    "CT": ("moderate", "Moderate protection (Ala/Val)", "Умеренная защита (Ala/Val)", 0.5),
                    "AC": ("moderate", "Moderate protection (Ala/Val)", "Умеренная защита (Ala/Val)", 0.5),
                    "TT": ("risk", "Lower mitochondrial protection (Val/Val)", "Сниженная митохондриальная защита (Val/Val)", -0.5),
                    "AA": ("risk", "Lower mitochondrial protection (Val/Val)", "Сниженная митохондриальная защита (Val/Val)", -0.5),
                    "AT": ("moderate", "Moderate protection (Ala/Val)", "Умеренная защита (Ala/Val)", 0.5),
                }
            },
            "rs1050450": {
                "gene": "GPX1",
                "description": "Glutathione peroxidase 1 - selenium-dependent antioxidant",
                "description_ru": "Глутатионпероксидаза 1 - селен-зависимый антиоксидант",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Pro/Pro - optimal GPX1 activity", "Pro/Pro - оптимальная активность GPX1", 0.8),
                    "CT": ("moderate", "Pro/Leu - moderate activity", "Pro/Leu - умеренная активность", 0.4),
                    "TT": ("risk", "Leu/Leu - reduced GPX1 activity", "Leu/Leu - сниженная активность GPX1", -0.4),
                }
            },
            "rs1001179": {
                "gene": "CAT",
                "description": "Catalase - hydrogen peroxide neutralization",
                "description_ru": "Каталаза - нейтрализация перекиси водорода",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Higher catalase activity", "Повышенная активность каталазы", 0.6),
                    "CT": ("moderate", "Average catalase", "Средняя активность каталазы", 0.3),
                    "TT": ("baseline", "Standard catalase", "Стандартная активность", 0),
                }
            },
        }
    },

    "inflammation": {
        "name": "Inflammaging Markers",
        "name_ru": "Маркеры воспалительного старения",
        "snps": {
            "rs1800795": {
                "gene": "IL-6",
                "description": "Interleukin-6 - pro-inflammatory cytokine",
                "description_ru": "Интерлейкин-6 - провоспалительный цитокин",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Lower IL-6 production - reduced inflammaging", "Низкая продукция IL-6 - меньше воспалительного старения", 1.2),
                    "CG": ("moderate", "Moderate IL-6 levels", "Умеренные уровни IL-6", 0.6),
                    "GG": ("risk", "Higher IL-6 - increased inflammaging", "Высокие уровни IL-6 - усиленное воспалительное старение", -0.8),
                }
            },
            "rs1800629": {
                "gene": "TNF-alpha",
                "description": "Tumor necrosis factor alpha - master inflammatory regulator",
                "description_ru": "Фактор некроза опухолей альфа - главный регулятор воспаления",
                "protective_allele": "G",
                "interpretation": {
                    "GG": ("beneficial", "Standard TNF-alpha production", "Стандартная продукция TNF-альфа", 0.5),
                    "GA": ("moderate", "Slightly elevated TNF-alpha", "Немного повышенный TNF-альфа", 0),
                    "AG": ("moderate", "Slightly elevated TNF-alpha", "Немного повышенный TNF-альфа", 0),
                    "AA": ("risk", "High TNF-alpha - chronic inflammation risk", "Высокий TNF-альфа - риск хронического воспаления", -1.0),
                }
            },
            "rs1205": {
                "gene": "CRP",
                "description": "C-reactive protein - inflammatory biomarker",
                "description_ru": "С-реактивный белок - биомаркер воспаления",
                "protective_allele": "T",
                "interpretation": {
                    "TT": ("beneficial", "Lower baseline CRP", "Низкий базовый уровень СРБ", 0.8),
                    "CT": ("moderate", "Average CRP levels", "Средние уровни СРБ", 0.4),
                    "CC": ("baseline", "Standard CRP", "Стандартные уровни СРБ", 0),
                }
            },
        }
    },

    "dna_repair": {
        "name": "DNA Repair",
        "name_ru": "Репарация ДНК",
        "snps": {
            "rs25487": {
                "gene": "XRCC1",
                "description": "X-ray repair cross-complementing 1 - base excision repair",
                "description_ru": "XRCC1 - эксцизионная репарация оснований",
                "protective_allele": "G",
                "interpretation": {
                    "GG": ("beneficial", "Arg/Arg - optimal DNA repair", "Arg/Arg - оптимальная репарация ДНК", 0.8),
                    "CC": ("beneficial", "Arg/Arg - optimal DNA repair", "Arg/Arg - оптимальная репарация ДНК", 0.8),
                    "AG": ("moderate", "Arg/Gln - moderate repair capacity", "Arg/Gln - умеренная способность к репарации", 0.4),
                    "GA": ("moderate", "Arg/Gln - moderate repair capacity", "Arg/Gln - умеренная способность к репарации", 0.4),
                    "CT": ("moderate", "Arg/Gln - moderate repair capacity", "Arg/Gln - умеренная способность к репарации", 0.4),
                    "TC": ("moderate", "Arg/Gln - moderate repair capacity", "Arg/Gln - умеренная способность к репарации", 0.4),
                    "AA": ("risk", "Gln/Gln - reduced DNA repair", "Gln/Gln - сниженная репарация ДНК", -0.6),
                    "TT": ("risk", "Gln/Gln - reduced DNA repair", "Gln/Gln - сниженная репарация ДНК", -0.6),
                }
            },
            "rs1052133": {
                "gene": "OGG1",
                "description": "8-oxoguanine DNA glycosylase - oxidative DNA damage repair",
                "description_ru": "OGG1 - репарация окислительных повреждений ДНК",
                "protective_allele": "C",
                "interpretation": {
                    "CC": ("beneficial", "Ser/Ser - efficient repair of oxidative damage", "Ser/Ser - эффективная репарация окислительных повреждений", 0.7),
                    "CG": ("moderate", "Ser/Cys - moderate repair", "Ser/Cys - умеренная репарация", 0.3),
                    "GG": ("risk", "Cys/Cys - reduced oxidative damage repair", "Cys/Cys - сниженная репарация окислительных повреждений", -0.5),
                }
            },
        }
    },

    "sirtuin": {
        "name": "Sirtuin Pathway",
        "name_ru": "Сиртуиновый путь",
        "snps": {
            "rs3758391": {
                "gene": "SIRT1",
                "description": "Sirtuin 1 - NAD+-dependent deacetylase, caloric restriction mimic",
                "description_ru": "Сиртуин 1 - NAD+-зависимая деацетилаза, имитация калорийного ограничения",
                "protective_allele": "T",
                "interpretation": {
                    "TT": ("beneficial", "Enhanced SIRT1 activity - better aging response", "Повышенная активность SIRT1 - лучший ответ на старение", 1.5),
                    "CT": ("moderate", "Moderate SIRT1 activity", "Умеренная активность SIRT1", 0.7),
                    "CC": ("baseline", "Standard SIRT1", "Стандартная активность SIRT1", 0),
                }
            },
        }
    },

    "igf1_pathway": {
        "name": "IGF-1/Insulin Pathway",
        "name_ru": "Путь IGF-1/Инсулин",
        "snps": {
            "rs2229765": {
                "gene": "IGF1R",
                "description": "Insulin-like growth factor 1 receptor - growth/longevity trade-off",
                "description_ru": "Рецептор инсулиноподобного фактора роста 1 - компромисс рост/долголетие",
                "protective_allele": "A",
                "interpretation": {
                    "AA": ("beneficial", "Reduced IGF-1R signaling - longevity pattern", "Сниженная сигнализация IGF-1R - паттерн долголетия", 1.3),
                    "AG": ("moderate", "Moderate IGF-1R signaling", "Умеренная сигнализация IGF-1R", 0.6),
                    "GA": ("moderate", "Moderate IGF-1R signaling", "Умеренная сигнализация IGF-1R", 0.6),
                    "GG": ("baseline", "Standard IGF-1R signaling", "Стандартная сигнализация IGF-1R", 0),
                }
            },
        }
    },
}

# APOE Genotype determination table
# rs429358: T=ancestral, C=derived (ε4)
# rs7412: C=ancestral, T=derived (ε2)
APOE_GENOTYPES = {
    ("TT", "CC"): ("ε3/ε3", "Most common genotype - baseline risk", "Самый распространённый генотип - базовый риск"),
    ("TT", "CT"): ("ε2/ε3", "Protective - lower cardiovascular and Alzheimer's risk", "Защитный - сниженный риск сердечно-сосудистых и Альцгеймера"),
    ("TT", "TT"): ("ε2/ε2", "Highly protective - lowest Alzheimer's risk, but increased triglycerides", "Высоко защитный - минимальный риск Альцгеймера, но повышены триглицериды"),
    ("CT", "CC"): ("ε3/ε4", "One ε4 allele - moderately increased Alzheimer's risk", "Один аллель ε4 - умеренно повышенный риск Альцгеймера"),
    ("CT", "CT"): ("ε2/ε4", "Mixed - ε2 partially compensates for ε4", "Смешанный - ε2 частично компенсирует ε4"),
    ("CC", "CC"): ("ε4/ε4", "Two ε4 alleles - significantly increased Alzheimer's risk", "Два аллеля ε4 - значительно повышенный риск Альцгеймера"),
}

# Anti-aging recommendations based on genotype patterns
RECOMMENDATIONS = {
    "oxidative_stress": {
        "risk": [
            "Increase antioxidant intake: Vitamin C (500-1000mg), Vitamin E (400 IU), CoQ10 (100-200mg)",
            "Consider NAC (N-Acetyl Cysteine) 600-1200mg/day for glutathione support",
            "Eat colorful vegetables rich in polyphenols",
            "Avoid excessive iron supplementation",
            "Consider astaxanthin supplementation (4-12mg/day)",
        ],
        "risk_ru": [
            "Увеличьте потребление антиоксидантов: Витамин C (500-1000мг), Витамин E (400 МЕ), CoQ10 (100-200мг)",
            "Рассмотрите NAC (N-ацетилцистеин) 600-1200мг/день для поддержки глутатиона",
            "Ешьте разноцветные овощи, богатые полифенолами",
            "Избегайте избыточного приёма железа",
            "Рассмотрите астаксантин (4-12мг/день)",
        ],
    },
    "inflammation": {
        "risk": [
            "Follow anti-inflammatory diet (Mediterranean, low glycemic)",
            "Consider omega-3 fatty acids (EPA/DHA 2-4g/day)",
            "Curcumin with piperine (500-1000mg/day)",
            "Reduce refined carbohydrates and processed foods",
            "Regular moderate exercise (avoid overtraining)",
            "Optimize sleep (7-9 hours)",
        ],
        "risk_ru": [
            "Придерживайтесь противовоспалительной диеты (средиземноморская, низкогликемическая)",
            "Рассмотрите омега-3 жирные кислоты (EPA/DHA 2-4г/день)",
            "Куркумин с пиперином (500-1000мг/день)",
            "Сократите рафинированные углеводы и переработанные продукты",
            "Регулярные умеренные упражнения (избегайте перетренированности)",
            "Оптимизируйте сон (7-9 часов)",
        ],
    },
    "telomeres": {
        "risk": [
            "Prioritize stress management (meditation, yoga)",
            "Ensure adequate sleep for telomere maintenance",
            "Consider TA-65 or astragalus root extract",
            "Regular moderate exercise",
            "Avoid excessive endurance training",
            "Optimize vitamin D levels (40-60 ng/mL)",
        ],
        "risk_ru": [
            "Приоритизируйте управление стрессом (медитация, йога)",
            "Обеспечьте достаточный сон для поддержания теломер",
            "Рассмотрите TA-65 или экстракт астрагала",
            "Регулярные умеренные упражнения",
            "Избегайте чрезмерных тренировок на выносливость",
            "Оптимизируйте уровень витамина D (40-60 нг/мл)",
        ],
    },
    "dna_repair": {
        "risk": [
            "Ensure adequate B vitamins (especially B12, folate)",
            "Consider NMN or NR for NAD+ support (250-500mg/day)",
            "Minimize UV and radiation exposure",
            "Avoid DNA-damaging substances (tobacco, excessive alcohol)",
            "Consider sulforaphane from broccoli sprouts",
        ],
        "risk_ru": [
            "Обеспечьте достаточное потребление витаминов группы B (особенно B12, фолат)",
            "Рассмотрите NMN или NR для поддержки NAD+ (250-500мг/день)",
            "Минимизируйте воздействие УФ и радиации",
            "Избегайте веществ, повреждающих ДНК (табак, избыточный алкоголь)",
            "Рассмотрите сульфорафан из ростков брокколи",
        ],
    },
    "sirtuin": {
        "baseline": [
            "Practice intermittent fasting or time-restricted eating",
            "Consider resveratrol (250-500mg/day) or pterostilbene",
            "NMN or NR supplementation for NAD+ boost",
            "Regular exercise activates sirtuins",
            "Cold exposure (cold showers, cryotherapy)",
        ],
        "baseline_ru": [
            "Практикуйте интервальное голодание или ограниченное по времени питание",
            "Рассмотрите ресвератрол (250-500мг/день) или птеростильбен",
            "Добавки NMN или NR для повышения NAD+",
            "Регулярные упражнения активируют сиртуины",
            "Холодовое воздействие (холодный душ, криотерапия)",
        ],
    },
    "apoe_e4": {
        "risk": [
            "CRITICAL: Prioritize cardiovascular and brain health",
            "Follow strict Mediterranean or MIND diet",
            "Regular aerobic exercise (150+ min/week)",
            "Optimize sleep and treat sleep apnea",
            "Monitor and control blood pressure, glucose, cholesterol",
            "Consider DHA supplementation (1-2g/day)",
            "Avoid head injuries",
            "Engage in cognitive activities and social connections",
            "Regular cardiovascular screening",
        ],
        "risk_ru": [
            "ВАЖНО: Приоритет здоровью сердца и мозга",
            "Строго следуйте средиземноморской или MIND диете",
            "Регулярные аэробные упражнения (150+ мин/неделю)",
            "Оптимизируйте сон и лечите апноэ сна",
            "Контролируйте давление, глюкозу, холестерин",
            "Рассмотрите DHA (1-2г/день)",
            "Избегайте травм головы",
            "Занимайтесь когнитивной активностью и поддерживайте социальные связи",
            "Регулярный сердечно-сосудистый скрининг",
        ],
    },
    "general_longevity": [
        "Caloric restriction or intermittent fasting",
        "Regular exercise combining cardio and resistance training",
        "Optimize sleep quality and circadian rhythm",
        "Stress management and social connections",
        "Regular health screenings",
        "Consider metformin or rapamycin (consult physician)",
        "Maintain healthy body weight",
        "Avoid smoking and limit alcohol",
    ],
    "general_longevity_ru": [
        "Калорийное ограничение или интервальное голодание",
        "Регулярные упражнения, сочетающие кардио и силовые",
        "Оптимизируйте качество сна и циркадный ритм",
        "Управление стрессом и социальные связи",
        "Регулярные медицинские обследования",
        "Рассмотрите метформин или рапамицин (консультация с врачом)",
        "Поддерживайте здоровый вес",
        "Избегайте курения и ограничьте алкоголь",
    ],
}


def load_genome():
    """Load and parse 23andMe genome file"""
    genome = {}
    with open(GENOME_FILE, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rsid = parts[0]
                genotype = parts[3]
                genome[rsid] = genotype
    return genome


def determine_apoe_genotype(genome):
    """
    Determine APOE genotype (ε2/ε2, ε2/ε3, ε3/ε3, ε2/ε4, ε3/ε4, ε4/ε4)
    Based on rs429358 and rs7412
    """
    rs429358 = genome.get("rs429358", "")
    rs7412 = genome.get("rs7412", "")

    # Normalize genotypes (sort alleles)
    def normalize(gt):
        if len(gt) == 2:
            return "".join(sorted(gt))
        return gt

    rs429358 = normalize(rs429358)
    rs7412 = normalize(rs7412)

    # Look up in table
    key = (rs429358, rs7412)

    if key in APOE_GENOTYPES:
        return {
            "genotype": APOE_GENOTYPES[key][0],
            "interpretation": APOE_GENOTYPES[key][1],
            "interpretation_ru": APOE_GENOTYPES[key][2],
            "rs429358": rs429358,
            "rs7412": rs7412,
        }
    else:
        # Manual determination based on individual SNPs
        # rs429358: C=ε4, T=not ε4
        # rs7412: T=ε2, C=not ε2
        e4_count = rs429358.count('C')
        e2_count = rs7412.count('T')

        if e4_count == 0 and e2_count == 0:
            genotype = "ε3/ε3"
        elif e4_count == 0 and e2_count == 1:
            genotype = "ε2/ε3"
        elif e4_count == 0 and e2_count == 2:
            genotype = "ε2/ε2"
        elif e4_count == 1 and e2_count == 0:
            genotype = "ε3/ε4"
        elif e4_count == 1 and e2_count == 1:
            genotype = "ε2/ε4"
        elif e4_count == 2:
            genotype = "ε4/ε4"
        else:
            genotype = "Unknown"

        return {
            "genotype": genotype,
            "interpretation": f"Determined from rs429358={rs429358}, rs7412={rs7412}",
            "interpretation_ru": f"Определено из rs429358={rs429358}, rs7412={rs7412}",
            "rs429358": rs429358,
            "rs7412": rs7412,
        }


def analyze_longevity(genome):
    """Analyze all longevity-related SNPs"""
    results = {}

    for category, category_data in LONGEVITY_SNPS.items():
        category_results = {
            "name": category_data["name"],
            "name_ru": category_data["name_ru"],
            "snps": {},
            "summary": {
                "beneficial": 0,
                "moderate": 0,
                "baseline": 0,
                "risk": 0,
            }
        }

        for rsid, snp_data in category_data["snps"].items():
            genotype = genome.get(rsid, "")

            if not genotype or genotype == "--":
                interpretation = {
                    "status": "not_tested",
                    "description": "Not tested in this chip",
                    "description_ru": "Не тестировался на этом чипе",
                    "score": 0,
                }
            else:
                # Normalize genotype
                normalized = "".join(sorted(genotype)) if len(genotype) == 2 else genotype

                # Try both original and normalized
                if genotype in snp_data["interpretation"]:
                    interp = snp_data["interpretation"][genotype]
                elif normalized in snp_data["interpretation"]:
                    interp = snp_data["interpretation"][normalized]
                else:
                    interp = ("unknown", f"Genotype {genotype} not in database", f"Генотип {genotype} не в базе", 0)

                interpretation = {
                    "status": interp[0],
                    "description": interp[1],
                    "description_ru": interp[2],
                    "score": interp[3] if len(interp) > 3 else 0,
                }

                if interp[0] in category_results["summary"]:
                    category_results["summary"][interp[0]] += 1

            category_results["snps"][rsid] = {
                "gene": snp_data["gene"],
                "description": snp_data["description"],
                "description_ru": snp_data["description_ru"],
                "genotype": genotype,
                "interpretation": interpretation,
            }

        results[category] = category_results

    # Add APOE genotype determination
    results["apoe_determination"] = determine_apoe_genotype(genome)

    return results


def calculate_longevity_score(results):
    """Calculate aggregate longevity score based on all analyzed SNPs"""
    total_score = 0
    scores_by_category = {}

    for category, data in results.items():
        if category == "apoe_determination":
            # Add APOE-specific scoring
            apoe = data["genotype"]
            if "ε2/ε2" in apoe:
                apoe_score = 2.5
            elif "ε2/ε3" in apoe:
                apoe_score = 1.5
            elif "ε3/ε3" in apoe:
                apoe_score = 0
            elif "ε2/ε4" in apoe:
                apoe_score = -0.5
            elif "ε3/ε4" in apoe:
                apoe_score = -2.0
            elif "ε4/ε4" in apoe:
                apoe_score = -4.0
            else:
                apoe_score = 0
            scores_by_category["apoe"] = apoe_score
            total_score += apoe_score
            continue

        if "snps" not in data:
            continue

        category_score = 0
        snp_count = 0

        for rsid, snp_data in data["snps"].items():
            if snp_data["interpretation"]["status"] != "not_tested":
                category_score += snp_data["interpretation"]["score"]
                snp_count += 1

        scores_by_category[category] = category_score
        total_score += category_score

    # Calculate percentile (rough estimate based on normal distribution)
    # Assuming average score is 0, std dev ~5
    import math

    def normal_cdf(x, mean=0, std=5):
        return 0.5 * (1 + math.erf((x - mean) / (std * math.sqrt(2))))

    percentile = normal_cdf(total_score) * 100

    return {
        "total_score": total_score,
        "category_scores": scores_by_category,
        "percentile": percentile,
        "interpretation": get_score_interpretation(total_score, percentile),
    }


def get_score_interpretation(score, percentile):
    """Get interpretation of longevity score"""
    if score >= 8:
        return {
            "level": "exceptional",
            "description": "Exceptional longevity genetics - top tier",
            "description_ru": "Исключительная генетика долголетия - высший уровень",
        }
    elif score >= 4:
        return {
            "level": "favorable",
            "description": "Favorable longevity profile",
            "description_ru": "Благоприятный профиль долголетия",
        }
    elif score >= 0:
        return {
            "level": "average",
            "description": "Average longevity genetics",
            "description_ru": "Средняя генетика долголетия",
        }
    elif score >= -4:
        return {
            "level": "below_average",
            "description": "Below average - lifestyle interventions important",
            "description_ru": "Ниже среднего - важны изменения образа жизни",
        }
    else:
        return {
            "level": "challenging",
            "description": "Challenging profile - proactive health management essential",
            "description_ru": "Сложный профиль - необходимо активное управление здоровьем",
        }


def get_recommendations(results, score_data):
    """Generate personalized anti-aging recommendations"""
    recommendations = []

    # Check each category for risk factors
    for category, data in results.items():
        if category == "apoe_determination":
            # APOE-specific recommendations
            if "ε4" in data["genotype"]:
                recommendations.append({
                    "category": "APOE ε4 Carrier",
                    "category_ru": "Носитель APOE ε4",
                    "priority": "high",
                    "items": RECOMMENDATIONS["apoe_e4"]["risk"],
                    "items_ru": RECOMMENDATIONS["apoe_e4"]["risk_ru"],
                })
            continue

        if "snps" not in data:
            continue

        # Count risks in category
        risk_count = data["summary"].get("risk", 0)

        if risk_count > 0 and category in RECOMMENDATIONS:
            cat_recs = RECOMMENDATIONS[category]
            if "risk" in cat_recs:
                recommendations.append({
                    "category": data["name"],
                    "category_ru": data["name_ru"],
                    "priority": "high" if risk_count > 1 else "moderate",
                    "items": cat_recs["risk"],
                    "items_ru": cat_recs["risk_ru"],
                })
        elif category in RECOMMENDATIONS and "baseline" in RECOMMENDATIONS[category]:
            # Add baseline recommendations for categories without risk
            recommendations.append({
                "category": data["name"],
                "category_ru": data["name_ru"],
                "priority": "low",
                "items": RECOMMENDATIONS[category]["baseline"],
                "items_ru": RECOMMENDATIONS[category]["baseline_ru"],
            })

    # Always add general longevity recommendations
    recommendations.append({
        "category": "General Longevity",
        "category_ru": "Общее долголетие",
        "priority": "moderate",
        "items": RECOMMENDATIONS["general_longevity"],
        "items_ru": RECOMMENDATIONS["general_longevity_ru"],
    })

    # Sort by priority
    priority_order = {"high": 0, "moderate": 1, "low": 2}
    recommendations.sort(key=lambda x: priority_order.get(x["priority"], 3))

    return recommendations


def generate_report(results, score_data, recommendations):
    """Generate markdown report"""
    report = []
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report.append("# Longevity & Aging Genetics Report")
    report.append(f"# Отчёт по генетике долголетия и старения\n")
    report.append(f"**Generated / Сгенерировано:** {timestamp}\n")
    report.append("---\n")

    # Executive Summary
    report.append("## Executive Summary / Краткое резюме\n")
    report.append(f"**Longevity Score / Балл долголетия:** {score_data['total_score']:.1f}")
    report.append(f"**Percentile / Процентиль:** {score_data['percentile']:.0f}%")
    report.append(f"**Assessment / Оценка:** {score_data['interpretation']['description']}")
    report.append(f"**Оценка (RU):** {score_data['interpretation']['description_ru']}\n")

    # Score breakdown
    report.append("### Score by Category / Баллы по категориям\n")
    report.append("| Category / Категория | Score / Балл |")
    report.append("|----------------------|--------------|")
    for cat, cat_score in score_data["category_scores"].items():
        cat_name = cat.replace("_", " ").title()
        score_str = f"+{cat_score:.1f}" if cat_score > 0 else f"{cat_score:.1f}"
        report.append(f"| {cat_name} | {score_str} |")
    report.append("")

    # APOE Genotype
    report.append("## APOE Genotype / Генотип APOE\n")
    apoe = results["apoe_determination"]
    report.append(f"**Genotype / Генотип:** {apoe['genotype']}")
    report.append(f"**Status / Статус:** {apoe['interpretation']}")
    report.append(f"**Статус (RU):** {apoe['interpretation_ru']}")
    report.append(f"- rs429358: {apoe['rs429358']}")
    report.append(f"- rs7412: {apoe['rs7412']}\n")

    if "ε4" in apoe["genotype"]:
        report.append("⚠️ **WARNING / ВНИМАНИЕ:** APOE ε4 carrier detected. See recommendations section.\n")
        report.append("⚠️ **ВНИМАНИЕ:** Обнаружен носитель APOE ε4. См. раздел рекомендаций.\n")
    elif "ε2" in apoe["genotype"]:
        report.append("✅ **FAVORABLE / БЛАГОПРИЯТНО:** APOE ε2 allele detected - protective against Alzheimer's.\n")
        report.append("✅ **БЛАГОПРИЯТНО:** Обнаружен аллель APOE ε2 - защита от болезни Альцгеймера.\n")

    # Detailed Results by Category
    report.append("---\n")
    report.append("## Detailed Analysis / Детальный анализ\n")

    for category, data in results.items():
        if category == "apoe_determination":
            continue
        if "snps" not in data:
            continue

        report.append(f"### {data['name']} / {data['name_ru']}\n")

        # Summary counts
        summary = data["summary"]
        status_icons = {"beneficial": "✅", "moderate": "➖", "baseline": "⚪", "risk": "⚠️"}
        summary_parts = []
        for status, count in summary.items():
            if count > 0:
                icon = status_icons.get(status, "")
                summary_parts.append(f"{icon} {status.title()}: {count}")
        if summary_parts:
            report.append(f"**Summary:** {' | '.join(summary_parts)}\n")

        # SNP table
        report.append("| SNP | Gene / Ген | Genotype / Генотип | Status / Статус | Interpretation / Интерпретация |")
        report.append("|-----|------------|-------------------|-----------------|-------------------------------|")

        for rsid, snp_data in data["snps"].items():
            gene = snp_data["gene"]
            genotype = snp_data["genotype"] or "N/A"
            status = snp_data["interpretation"]["status"]
            desc_ru = snp_data["interpretation"]["description_ru"]

            status_icon = status_icons.get(status, "")
            report.append(f"| {rsid} | {gene} | {genotype} | {status_icon} {status} | {desc_ru} |")

        report.append("")

    # Recommendations
    report.append("---\n")
    report.append("## Anti-Aging Recommendations / Рекомендации по антистарению\n")

    priority_labels = {
        "high": "🔴 HIGH PRIORITY / ВЫСОКИЙ ПРИОРИТЕТ",
        "moderate": "🟡 MODERATE / УМЕРЕННЫЙ",
        "low": "🟢 GENERAL / ОБЩИЕ",
    }

    for rec in recommendations:
        priority = priority_labels.get(rec["priority"], rec["priority"])
        report.append(f"### {rec['category']} / {rec['category_ru']}")
        report.append(f"**Priority / Приоритет:** {priority}\n")

        for item, item_ru in zip(rec["items"], rec["items_ru"]):
            report.append(f"- {item}")
            report.append(f"  - *{item_ru}*")
        report.append("")

    # Disclaimer
    report.append("---\n")
    report.append("## Disclaimer / Отказ от ответственности\n")
    report.append("""
This report is for educational and informational purposes only. It is NOT medical advice.
Genetic factors are only part of the longevity equation. Lifestyle, environment, and other
factors play significant roles. Consult healthcare professionals before making health decisions.

Этот отчёт предназначен только для образовательных и информационных целей. Это НЕ медицинский совет.
Генетические факторы - лишь часть уравнения долголетия. Образ жизни, окружающая среда и другие
факторы играют значительную роль. Консультируйтесь с медицинскими специалистами перед принятием решений о здоровье.
""")

    return "\n".join(report)


def main():
    """Main execution function"""
    print("=" * 60)
    print("LONGEVITY & AGING GENETICS ANALYSIS")
    print("АНАЛИЗ ГЕНЕТИКИ ДОЛГОЛЕТИЯ И СТАРЕНИЯ")
    print("=" * 60)
    print()

    # Load genome
    print("Loading genome data...")
    genome = load_genome()
    print(f"Loaded {len(genome):,} SNPs from genome file.\n")

    # Analyze longevity SNPs
    print("Analyzing longevity markers...")
    results = analyze_longevity(genome)

    # Calculate longevity score
    print("Calculating longevity score...")
    score_data = calculate_longevity_score(results)

    # Get recommendations
    print("Generating personalized recommendations...")
    recommendations = get_recommendations(results, score_data)

    # Generate report
    print("Generating report...")
    report = generate_report(results, score_data, recommendations)

    # Save report
    report_dir = f"{REPORTS_PATH}/longevity"
    os.makedirs(report_dir, exist_ok=True)
    report_file = f"{report_dir}/report.md"

    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)

    print(f"\nReport saved to: {report_file}")

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY / РЕЗЮМЕ")
    print("=" * 60)
    print(f"\nLongevity Score: {score_data['total_score']:.1f}")
    print(f"Percentile: {score_data['percentile']:.0f}%")
    print(f"Assessment: {score_data['interpretation']['description']}")
    print(f"Оценка: {score_data['interpretation']['description_ru']}")

    print(f"\nAPOE Genotype: {results['apoe_determination']['genotype']}")
    print(f"Status: {results['apoe_determination']['interpretation']}")

    print("\n" + "-" * 60)
    print("Category Scores:")
    for cat, cat_score in score_data["category_scores"].items():
        score_str = f"+{cat_score:.1f}" if cat_score > 0 else f"{cat_score:.1f}"
        print(f"  {cat.replace('_', ' ').title()}: {score_str}")

    print("\n" + "=" * 60)
    print("Analysis complete. See report for detailed recommendations.")
    print("=" * 60)


if __name__ == "__main__":
    main()
