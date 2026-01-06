#!/usr/bin/env python3
"""
Skin Analysis Script
Analyzes skin-related genetic markers from 23andMe data
"""

import os
from datetime import datetime
from collections import defaultdict

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Skin markers
# =============================================================================

SKIN_SNPS = {
    "collagen": {
        "name": "Коллаген и структура кожи",
        "snps": {
            "rs1800012": {
                "gene": "COL1A1",
                "description": "Структура коллагена I типа, склонность к морщинам",
                "interpretation": {
                    "GG": ("normal", "Нормальная структура коллагена"),
                    "GT": ("moderate", "Умеренно сниженная плотность коллагена"),
                    "TT": ("impaired", "Сниженная плотность коллагена, склонность к морщинам"),
                    "CC": ("normal", "Нормальная структура коллагена"),
                    "CT": ("moderate", "Умеренно сниженная плотность коллагена"),
                }
            },
        }
    },

    "mmp1": {
        "name": "Деградация коллагена (MMP1)",
        "snps": {
            "rs1799750": {
                "gene": "MMP1",
                "description": "Коллагеназа - расщепление коллагена",
                "interpretation": {
                    # 2G/2G variant - insertion/deletion polymorphism
                    # In 23andMe: usually reported as G/G, GG, or deletion status
                    "GG": ("fast_aging", "2G/2G - Быстрое разрушение коллагена, ранние морщины"),
                    "DG": ("moderate", "1G/2G - Умеренная деградация коллагена"),
                    "GD": ("moderate", "1G/2G - Умеренная деградация коллагена"),
                    "DD": ("normal", "1G/1G - Лучшее сохранение коллагена"),
                    # Alternative representations
                    "II": ("fast_aging", "2G/2G - Быстрое разрушение коллагена"),
                    "DI": ("moderate", "1G/2G - Умеренная деградация"),
                    "ID": ("moderate", "1G/2G - Умеренная деградация"),
                    # Some chips report as T/C
                    "TT": ("normal", "1G/1G - Лучшее сохранение коллагена"),
                    "CT": ("moderate", "1G/2G - Умеренная деградация"),
                    "TC": ("moderate", "1G/2G - Умеренная деградация"),
                    "CC": ("fast_aging", "2G/2G - Быстрое разрушение коллагена"),
                }
            },
        }
    },

    "antioxidants": {
        "name": "Антиоксидантная защита",
        "snps": {
            "rs4880": {
                "gene": "SOD2",
                "description": "Супероксиддисмутаза - защита от окислительного стресса",
                "interpretation": {
                    "CC": ("good", "Высокая антиоксидантная защита (Ala/Ala)"),
                    "CT": ("moderate", "Умеренная защита (Ala/Val)"),
                    "TC": ("moderate", "Умеренная защита (Ala/Val)"),
                    "TT": ("low", "Сниженная антиоксидантная защита (Val/Val)"),
                    "AA": ("good", "Высокая антиоксидантная защита"),
                    "AG": ("moderate", "Умеренная защита"),
                    "GA": ("moderate", "Умеренная защита"),
                    "GG": ("low", "Сниженная антиоксидантная защита"),
                }
            },
            "rs1695": {
                "gene": "GSTP1",
                "description": "Глутатион-S-трансфераза - детоксикация",
                "interpretation": {
                    "AA": ("good", "Высокая детоксикационная способность (Ile/Ile)"),
                    "AG": ("moderate", "Умеренная детоксикация (Ile/Val)"),
                    "GA": ("moderate", "Умеренная детоксикация (Ile/Val)"),
                    "GG": ("low", "Сниженная детоксикация (Val/Val)"),
                }
            },
        }
    },

    "uv_sensitivity": {
        "name": "Чувствительность к УФ-излучению",
        "snps": {
            "rs1805007": {
                "gene": "MC1R (R151C)",
                "description": "Рецептор меланокортина - рыжие волосы, веснушки",
                "interpretation": {
                    "CC": ("normal", "Обычная чувствительность к УФ"),
                    "CT": ("sensitive", "Носитель - повышенная чувствительность к солнцу"),
                    "TC": ("sensitive", "Носитель - повышенная чувствительность к солнцу"),
                    "TT": ("high_risk", "Высокая чувствительность к УФ, риск ожогов"),
                }
            },
            "rs1805008": {
                "gene": "MC1R (R160W)",
                "description": "Рецептор меланокортина - фоточувствительность",
                "interpretation": {
                    "CC": ("normal", "Обычная чувствительность к УФ"),
                    "CT": ("sensitive", "Носитель - повышенная чувствительность"),
                    "TC": ("sensitive", "Носитель - повышенная чувствительность"),
                    "TT": ("high_risk", "Высокая чувствительность, легко обгорает"),
                }
            },
            "rs12913832": {
                "gene": "HERC2/OCA2",
                "description": "Цвет глаз, пигментация кожи",
                "interpretation": {
                    "AA": ("light", "Светлые глаза, светлая кожа - выше чувствительность к УФ"),
                    "AG": ("medium", "Средняя пигментация"),
                    "GA": ("medium", "Средняя пигментация"),
                    "GG": ("dark", "Тёмные глаза, лучшая защита от УФ"),
                }
            },
        }
    },

    "photoaging": {
        "name": "Фотостарение",
        "snps": {
            "rs1805005": {
                "gene": "MC1R (V60L)",
                "description": "Рецептор меланокортина - фотостарение",
                "interpretation": {
                    "GG": ("normal", "Обычный риск фотостарения"),
                    "GT": ("elevated", "Повышенный риск фотостарения"),
                    "TG": ("elevated", "Повышенный риск фотостарения"),
                    "TT": ("high", "Высокий риск фотостарения"),
                    "AA": ("normal", "Обычный риск фотостарения"),
                    "AG": ("elevated", "Повышенный риск фотостарения"),
                    "GA": ("elevated", "Повышенный риск фотостарения"),
                }
            },
            "rs1805009": {
                "gene": "MC1R (D294H)",
                "description": "Рецептор меланокортина - пигментные пятна",
                "interpretation": {
                    "GG": ("normal", "Обычный риск пигментации"),
                    "GA": ("elevated", "Повышенный риск пигментных пятен"),
                    "AG": ("elevated", "Повышенный риск пигментных пятен"),
                    "AA": ("high", "Высокий риск гиперпигментации"),
                    "CC": ("normal", "Обычный риск пигментации"),
                    "CA": ("elevated", "Повышенный риск пигментных пятен"),
                    "AC": ("elevated", "Повышенный риск пигментных пятен"),
                }
            },
        }
    },

    "acne": {
        "name": "Акне и воспаления кожи",
        "snps": {
            "rs4133274": {
                "gene": "DDB2",
                "description": "Репарация ДНК, склонность к акне",
                "interpretation": {
                    "TT": ("normal", "Обычный риск акне"),
                    "TG": ("elevated", "Умеренно повышенный риск акне"),
                    "GT": ("elevated", "Умеренно повышенный риск акне"),
                    "GG": ("high", "Повышенный риск акне"),
                    "CC": ("normal", "Обычный риск акне"),
                    "CA": ("elevated", "Умеренно повышенный риск"),
                    "AC": ("elevated", "Умеренно повышенный риск"),
                    "AA": ("high", "Повышенный риск акне"),
                }
            },
            "rs1800629": {
                "gene": "TNF-alpha",
                "description": "Фактор некроза опухоли - воспаление",
                "interpretation": {
                    "GG": ("normal", "Нормальный уровень воспаления"),
                    "GA": ("elevated", "Повышенная склонность к воспалению"),
                    "AG": ("elevated", "Повышенная склонность к воспалению"),
                    "AA": ("high", "Высокая склонность к воспалению кожи"),
                }
            },
        }
    },

    "psoriasis": {
        "name": "Псориаз",
        "snps": {
            "rs10484554": {
                "gene": "HLA-C",
                "description": "Главный генетический фактор риска псориаза",
                "interpretation": {
                    "CC": ("normal", "Низкий генетический риск псориаза"),
                    "CT": ("elevated", "Умеренный генетический риск псориаза"),
                    "TC": ("elevated", "Умеренный генетический риск псориаза"),
                    "TT": ("high", "Высокий генетический риск псориаза"),
                }
            },
        }
    },

    "eczema": {
        "name": "Экзема (атопический дерматит)",
        "snps": {
            "rs61816761": {
                "gene": "FLG (Filaggrin)",
                "description": "Филаггрин - барьерная функция кожи",
                "interpretation": {
                    "GG": ("normal", "Нормальный кожный барьер"),
                    "GA": ("impaired", "Нарушение кожного барьера, риск экземы"),
                    "AG": ("impaired", "Нарушение кожного барьера, риск экземы"),
                    "AA": ("high_risk", "Высокий риск атопического дерматита"),
                    "CC": ("normal", "Нормальный кожный барьер"),
                    "CT": ("impaired", "Нарушение кожного барьера"),
                    "TC": ("impaired", "Нарушение кожного барьера"),
                    "TT": ("high_risk", "Высокий риск экземы"),
                }
            },
        }
    },

    "wound_healing": {
        "name": "Заживление ран",
        "snps": {
            "rs1800629_wound": {
                "gene": "TNF-alpha",
                "snp_id_actual": "rs1800629",
                "description": "Скорость заживления, рубцевание",
                "interpretation": {
                    "GG": ("normal", "Нормальное заживление"),
                    "GA": ("slow", "Замедленное заживление, риск келоидов"),
                    "AG": ("slow", "Замедленное заживление, риск келоидов"),
                    "AA": ("impaired", "Склонность к плохому заживлению"),
                }
            },
            "rs1800795": {
                "gene": "IL-6",
                "description": "Интерлейкин-6 - воспаление и заживление",
                "interpretation": {
                    "GG": ("normal", "Нормальный воспалительный ответ"),
                    "GC": ("elevated", "Повышенное воспаление, медленнее заживление"),
                    "CG": ("elevated", "Повышенное воспаление, медленнее заживление"),
                    "CC": ("high", "Высокое воспаление, замедленное заживление"),
                }
            },
        }
    },

    "cellulite": {
        "name": "Целлюлит",
        "snps": {
            "rs1799750_cellulite": {
                "gene": "MMP1",
                "snp_id_actual": "rs1799750",
                "description": "Структура соединительной ткани",
                "interpretation": {
                    "GG": ("high_risk", "2G/2G - Повышенный риск целлюлита"),
                    "DG": ("moderate", "1G/2G - Умеренный риск"),
                    "GD": ("moderate", "1G/2G - Умеренный риск"),
                    "DD": ("normal", "1G/1G - Сниженный риск целлюлита"),
                    "TT": ("normal", "Сниженный риск целлюлита"),
                    "CT": ("moderate", "Умеренный риск"),
                    "TC": ("moderate", "Умеренный риск"),
                    "CC": ("high_risk", "Повышенный риск целлюлита"),
                }
            },
            "rs1800012_cellulite": {
                "gene": "COL1A1",
                "snp_id_actual": "rs1800012",
                "description": "Плотность коллагена, эластичность кожи",
                "interpretation": {
                    "GG": ("normal", "Нормальная плотность коллагена"),
                    "GT": ("moderate", "Умеренно снижена плотность"),
                    "TT": ("high_risk", "Сниженная плотность, риск целлюлита"),
                    "CC": ("normal", "Нормальная плотность коллагена"),
                    "CT": ("moderate", "Умеренно снижена плотность"),
                }
            },
        }
    },

    "elasticity": {
        "name": "Эластичность кожи",
        "snps": {
            "rs7539120": {
                "gene": "ELN (Elastin)",
                "description": "Эластин - упругость кожи",
                "interpretation": {
                    "CC": ("good", "Хорошая эластичность кожи"),
                    "CT": ("moderate", "Умеренная эластичность"),
                    "TC": ("moderate", "Умеренная эластичность"),
                    "TT": ("reduced", "Сниженная эластичность, ранние морщины"),
                    "GG": ("good", "Хорошая эластичность кожи"),
                    "GA": ("moderate", "Умеренная эластичность"),
                    "AG": ("moderate", "Умеренная эластичность"),
                    "AA": ("reduced", "Сниженная эластичность"),
                }
            },
        }
    },

    "hydration": {
        "name": "Увлажнённость кожи",
        "snps": {
            "rs12212041": {
                "gene": "AQP3",
                "description": "Аквапорин-3 - водный баланс кожи",
                "interpretation": {
                    "TT": ("good", "Хорошее удержание влаги"),
                    "TC": ("moderate", "Умеренное увлажнение"),
                    "CT": ("moderate", "Умеренное увлажнение"),
                    "CC": ("dry", "Склонность к сухой коже"),
                    "AA": ("good", "Хорошее удержание влаги"),
                    "AG": ("moderate", "Умеренное увлажнение"),
                    "GA": ("moderate", "Умеренное увлажнение"),
                    "GG": ("dry", "Склонность к сухой коже"),
                }
            },
        }
    },

    "glycation": {
        "name": "Гликация (повреждение сахарами)",
        "snps": {
            "rs2070600": {
                "gene": "AGER (RAGE)",
                "description": "Рецептор конечных продуктов гликирования",
                "interpretation": {
                    "CC": ("protective", "Защита от гликационного старения"),
                    "CT": ("moderate", "Умеренный риск гликации"),
                    "TC": ("moderate", "Умеренный риск гликации"),
                    "TT": ("elevated", "Повышенный риск гликационного старения"),
                    "GG": ("protective", "Защита от гликации"),
                    "GA": ("moderate", "Умеренный риск"),
                    "AG": ("moderate", "Умеренный риск"),
                    "AA": ("elevated", "Повышенный риск"),
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


def analyze_skin(genome):
    """Analyze skin markers"""
    results = {}

    for category, cat_info in SKIN_SNPS.items():
        cat_results = []
        for snp_id, snp_info in cat_info['snps'].items():
            # Handle duplicate SNPs with different interpretations
            actual_snp_id = snp_info.get('snp_id_actual', snp_id)

            result = {
                'snp_id': actual_snp_id,
                'gene': snp_info['gene'],
                'description': snp_info['description'],
                'found': False,
                'genotype': None,
                'status': None,
                'interpretation': None
            }

            if actual_snp_id in genome:
                result['found'] = True
                genotype = genome[actual_snp_id]['genotype']
                result['genotype'] = genotype

                interp = snp_info.get('interpretation', {})
                for gt in [genotype, genotype[::-1] if len(genotype) == 2 else genotype]:
                    if gt in interp:
                        result['status'], result['interpretation'] = interp[gt]
                        break

            cat_results.append(result)
        results[category] = cat_results

    return results


def determine_skin_profile(results):
    """Determine overall skin profile based on genetic results"""

    scores = {
        'aging_prone': 0,
        'sensitive': 0,
        'resilient': 0,
        'inflammation_prone': 0,
        'photoaging_risk': 0,
    }

    profile_factors = []

    # Analyze collagen and aging markers
    for category in ['collagen', 'mmp1', 'elasticity']:
        if category in results:
            for r in results[category]:
                if r['status'] in ['impaired', 'fast_aging', 'reduced']:
                    scores['aging_prone'] += 2
                    profile_factors.append(f"Быстрое старение ({r['gene']})")
                elif r['status'] in ['moderate']:
                    scores['aging_prone'] += 1
                elif r['status'] in ['normal', 'good']:
                    scores['resilient'] += 1

    # Analyze UV sensitivity and photoaging
    for category in ['uv_sensitivity', 'photoaging']:
        if category in results:
            for r in results[category]:
                if r['status'] in ['high_risk', 'high', 'sensitive', 'light']:
                    scores['sensitive'] += 2
                    scores['photoaging_risk'] += 2
                    profile_factors.append(f"УФ-чувствительность ({r['gene']})")
                elif r['status'] in ['elevated', 'medium']:
                    scores['sensitive'] += 1
                    scores['photoaging_risk'] += 1
                elif r['status'] in ['normal', 'dark']:
                    scores['resilient'] += 1

    # Analyze inflammation and skin conditions
    for category in ['acne', 'psoriasis', 'eczema', 'wound_healing']:
        if category in results:
            for r in results[category]:
                if r['status'] in ['high', 'high_risk', 'impaired']:
                    scores['inflammation_prone'] += 2
                    profile_factors.append(f"Склонность к воспалению ({r['gene']})")
                elif r['status'] in ['elevated', 'slow']:
                    scores['inflammation_prone'] += 1
                elif r['status'] in ['normal']:
                    scores['resilient'] += 1

    # Analyze antioxidant protection
    if 'antioxidants' in results:
        for r in results['antioxidants']:
            if r['status'] in ['low']:
                scores['aging_prone'] += 1
                scores['sensitive'] += 1
            elif r['status'] in ['good']:
                scores['resilient'] += 2

    # Determine primary profile
    profiles = []

    if scores['aging_prone'] >= 3:
        profiles.append("Склонность к раннему старению")
    if scores['sensitive'] >= 3 or scores['photoaging_risk'] >= 3:
        profiles.append("Чувствительная кожа")
    if scores['inflammation_prone'] >= 3:
        profiles.append("Склонность к воспалениям")
    if scores['resilient'] >= 5 and not profiles:
        profiles.append("Устойчивая кожа")

    if not profiles:
        profiles.append("Средняя устойчивость")

    return {
        'profiles': profiles,
        'scores': scores,
        'factors': profile_factors
    }


def generate_report(results, profile):
    """Generate skin analysis report"""
    report = []
    report.append("# Генетический анализ кожи")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("\n---\n")

    # Skin profile summary
    report.append("## Ваш генетический профиль кожи\n")

    profile_emoji = {
        'Склонность к раннему старению': '⏳',
        'Чувствительная кожа': '🌡️',
        'Склонность к воспалениям': '🔥',
        'Устойчивая кожа': '💪',
        'Средняя устойчивость': '⚖️',
    }

    for p in profile['profiles']:
        emoji = profile_emoji.get(p, '•')
        report.append(f"### {emoji} {p}\n")

    if profile['factors']:
        report.append("**Ключевые генетические факторы:**\n")
        for f in profile['factors'][:5]:  # Top 5 factors
            report.append(f"- {f}")
        report.append("")

    report.append("---\n")

    # Summary of key findings
    report.append("## Ключевые находки\n")

    warnings = []
    positive = []

    for category, cat_results in results.items():
        for r in cat_results:
            if r['found'] and r['status']:
                if r['status'] in ['impaired', 'fast_aging', 'high_risk', 'high', 'low', 'reduced', 'dry', 'elevated']:
                    warnings.append(f"⚠️ **{r['gene']}**: {r['interpretation']}")
                elif r['status'] in ['good', 'protective', 'normal', 'dark']:
                    if r['status'] in ['good', 'protective']:
                        positive.append(f"✅ **{r['gene']}**: {r['interpretation']}")

    if warnings:
        report.append("### Требуют внимания\n")
        for w in warnings:
            report.append(f"- {w}")
        report.append("")

    if positive:
        report.append("### Генетические преимущества\n")
        for p in positive:
            report.append(f"- {p}")
        report.append("")

    report.append("---\n")

    # Detailed results by category
    category_order = [
        'collagen', 'mmp1', 'elasticity', 'antioxidants',
        'uv_sensitivity', 'photoaging', 'hydration', 'glycation',
        'acne', 'psoriasis', 'eczema', 'wound_healing', 'cellulite'
    ]

    for category in category_order:
        if category not in results:
            continue
        cat_results = results[category]
        cat_name = SKIN_SNPS[category]['name']
        report.append(f"## {cat_name}\n")
        report.append("| SNP | Ген | Генотип | Статус | Интерпретация |")
        report.append("|-----|-----|---------|--------|---------------|")

        for r in cat_results:
            if r['found']:
                status_emoji = {
                    'normal': '✅',
                    'good': '✅',
                    'protective': '✅',
                    'dark': '✅',
                    'moderate': '🟡',
                    'medium': '🟡',
                    'slow': '🟡',
                    'elevated': '🟠',
                    'impaired': '🔴',
                    'fast_aging': '🔴',
                    'high_risk': '🔴',
                    'high': '🔴',
                    'low': '🔴',
                    'reduced': '🔴',
                    'dry': '🔴',
                    'sensitive': '🟠',
                    'light': 'ℹ️',
                }.get(r['status'], '•')
                interp = r['interpretation'] or 'Нет данных'
                status = r['status'] or 'н/д'
                report.append(f"| {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {status_emoji} {status} | {interp} |")
            else:
                report.append(f"| {r['snp_id']} | {r['gene']} | - | - | Не найден |")
        report.append("")

    # Recommendations section
    report.append("---\n")
    report.append("## Персонализированные рекомендации по уходу за кожей\n")

    # Generate recommendations based on profile
    recommendations = []

    # Aging-prone skin recommendations
    if profile['scores']['aging_prone'] >= 2:
        recommendations.append("""### Антивозрастной уход

- **Ретиноиды**: Начните с 0.025-0.05% ретинола 2-3 раза в неделю
- **Пептиды**: Ищите Matrixyl, Argireline для стимуляции коллагена
- **Витамин C**: 10-20% L-аскорбиновая кислота утром
- **Коллагеновые добавки**: 5-10г гидролизованного коллагена в день
- **Профилактика**: Избегайте курения и сахара (гликация)""")

    # UV-sensitive skin recommendations
    if profile['scores']['sensitive'] >= 2 or profile['scores']['photoaging_risk'] >= 2:
        recommendations.append("""### Защита от солнца

- **SPF 50+**: Ежедневно, даже в пасмурную погоду
- **Физические фильтры**: Оксид цинка, диоксид титана (менее раздражающие)
- **Обновление**: Каждые 2 часа при нахождении на солнце
- **Одежда**: UPF-одежда, широкополые шляпы, солнцезащитные очки
- **Время**: Избегайте прямого солнца с 10:00 до 16:00
- **Антиоксиданты**: Витамин C, E, ниацинамид для защиты от фотоповреждений""")

    # Inflammation-prone skin recommendations
    if profile['scores']['inflammation_prone'] >= 2:
        recommendations.append("""### Успокаивающий уход

- **Барьерные средства**: Церамиды, холестерол, жирные кислоты
- **Успокаивающие ингредиенты**: Ниацинамид (B3), центелла азиатская, аллантоин
- **Избегайте**: Агрессивных ПАВ, спирта, отдушек
- **Омега-3**: 2-3г EPA/DHA ежедневно (противовоспалительный эффект)
- **Пробиотики**: Для здоровья микробиома кожи
- **Диета**: Ограничьте молочные продукты и сахар при акне""")

    # Antioxidant support
    low_antioxidants = False
    if 'antioxidants' in results:
        for r in results['antioxidants']:
            if r['status'] in ['low', 'moderate']:
                low_antioxidants = True
                break

    if low_antioxidants:
        recommendations.append("""### Антиоксидантная поддержка

- **Витамин C**: 15-20% сыворотка утром
- **Витамин E**: В комбинации с витамином C
- **Ниацинамид**: 5-10% для защиты и восстановления
- **Ресвератрол**: Антиоксидант из винограда
- **Добавки**: Астаксантин 4-12мг, CoQ10 100-200мг
- **Питание**: Больше ягод, тёмной зелени, орехов""")

    # Hydration support
    if 'hydration' in results:
        for r in results['hydration']:
            if r['status'] in ['dry', 'moderate']:
                recommendations.append("""### Увлажнение

- **Гиалуроновая кислота**: Разные молекулярные веса
- **Церамиды**: Восстановление барьера
- **Глицерин, мочевина**: Увлажнители
- **Окклюзивы**: Сквалан, масло жожоба на ночь
- **Увлажнитель воздуха**: Особенно зимой
- **Вода**: 2-2.5л в день""")
                break

    # Glycation protection
    if 'glycation' in results:
        for r in results['glycation']:
            if r['status'] in ['elevated', 'moderate']:
                recommendations.append("""### Защита от гликации

- **Диета**: Ограничьте сахар и быстрые углеводы
- **Карнозин**: Добавка 500-1000мг в день
- **Альфа-липоевая кислота**: 100-300мг в день
- **Избегайте**: Жареной пищи с корочкой (AGEs)
- **Зелёный чай**: Полифенолы против гликации""")
                break

    if not recommendations:
        recommendations.append("""### Базовый уход

Ваши гены показывают хорошую устойчивость кожи. Придерживайтесь базового ухода:
- Мягкое очищение
- Увлажнение
- SPF 30+ ежедневно
- Антиоксидантная сыворотка""")

    for rec in recommendations:
        report.append(rec)
        report.append("")

    report.append("---\n")
    report.append("## Важно\n")
    report.append("- Генетика — это предрасположенность, не диагноз")
    report.append("- Образ жизни, диета и уход влияют на состояние кожи не меньше генов")
    report.append("- При кожных заболеваниях обратитесь к дерматологу")
    report.append("- Начинайте новые активные ингредиенты постепенно")

    return '\n'.join(report)


def main():
    print("=" * 60)
    print("ГЕНЕТИЧЕСКИЙ АНАЛИЗ КОЖИ")
    print("=" * 60)

    print("\n[1/4] Загрузка генома...")
    genome = load_genome()
    print(f"      Загружено {len(genome)} SNP")

    print("\n[2/4] Анализ маркеров кожи...")
    results = analyze_skin(genome)

    total = sum(len(r) for r in results.values())
    found = sum(sum(1 for x in r if x['found']) for r in results.values())
    print(f"      Найдено: {found}/{total} маркеров")

    print("\n[3/4] Определение профиля кожи...")
    profile = determine_skin_profile(results)
    print(f"      Профиль: {', '.join(profile['profiles'])}")

    print("\n[4/4] Генерация отчёта...")
    report = generate_report(results, profile)

    # Ensure directory exists
    report_dir = f"{REPORTS_PATH}/skin"
    os.makedirs(report_dir, exist_ok=True)

    report_path = f"{report_dir}/report.md"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"      -> {report_path}")

    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)


if __name__ == "__main__":
    main()
