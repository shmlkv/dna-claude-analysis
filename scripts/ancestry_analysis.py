#!/usr/bin/env python3
"""
Ancestry Analysis Script
Analyzes ancestry-related genetic markers from 23andMe data
"""

import os
from datetime import datetime

# Paths
BASE_PATH = "/Users/sh/Library/Mobile Documents/com~apple~CloudDocs/dna"
GENOME_FILE = f"{BASE_PATH}/data/genome_Andre_Sh_v5_Full_20260106100611.txt"
REPORTS_PATH = f"{BASE_PATH}/reports"

# =============================================================================
# SNP DATABASE - Ancestry markers
# =============================================================================

ANCESTRY_SNPS = {
    "pigmentation": {
        "name": "Пигментация (происхождение)",
        "snps": {
            "rs1426654": {
                "gene": "SLC24A5",
                "description": "Светлая кожа (Европа)",
                "interpretation": {
                    "AA": ("european", "Европейский вариант - светлая кожа"),
                    "AG": ("mixed", "Смешанный"),
                    "GG": ("ancestral", "Предковый вариант - тёмная кожа"),
                }
            },
            "rs16891982": {
                "gene": "SLC45A2",
                "description": "Светлая кожа (Европа)",
                "interpretation": {
                    "GG": ("european", "Европейский вариант - светлая кожа"),
                    "CG": ("mixed", "Смешанный"),
                    "CC": ("ancestral", "Предковый вариант"),
                }
            },
            "rs12913832": {
                "gene": "HERC2/OCA2",
                "description": "Цвет глаз, Северная Европа",
                "interpretation": {
                    "AA": ("northern_european", "Голубые глаза - Северная Европа"),
                    "AG": ("mixed", "Зелёные/ореховые глаза"),
                    "GG": ("other", "Карие глаза"),
                }
            },
            "rs1800407": {
                "gene": "OCA2",
                "description": "Светлые глаза (Европа)",
                "interpretation": {
                    "TT": ("european", "Светлые глаза"),
                    "GT": ("mixed", "Смешанный"),
                    "GG": ("other", "Тёмные глаза"),
                }
            },
        }
    },

    "east_asian": {
        "name": "Восточноазиатские маркеры",
        "snps": {
            "rs3827760": {
                "gene": "EDAR",
                "description": "Толстые волосы, лопатообразные резцы",
                "interpretation": {
                    "CC": ("east_asian", "Восточноазиатский вариант - толстые волосы"),
                    "CT": ("mixed", "Смешанный"),
                    "TT": ("other", "Не-азиатский вариант"),
                }
            },
            "rs17822931": {
                "gene": "ABCC11",
                "description": "Тип ушной серы",
                "interpretation": {
                    "CC": ("east_asian", "Сухая ушная сера (Восточная Азия)"),
                    "CT": ("mixed", "Смешанный"),
                    "TT": ("other", "Влажная ушная сера (Европа/Африка)"),
                }
            },
        }
    },

    "african": {
        "name": "Африканские маркеры",
        "snps": {
            "rs2814778": {
                "gene": "DARC (Duffy)",
                "description": "Duffy null - резистентность к малярии",
                "interpretation": {
                    "CC": ("african", "Duffy null - Африканское происхождение, защита от P.vivax"),
                    "CT": ("mixed", "Носитель"),
                    "TT": ("other", "Duffy positive"),
                }
            },
        }
    },

    "lactase": {
        "name": "Персистенция лактазы (география)",
        "snps": {
            "rs4988235": {
                "gene": "LCT (MCM6)",
                "description": "Переносимость лактозы",
                "interpretation": {
                    "TT": ("european_pastoral", "Европейский пастушеский вариант - переносит лактозу"),
                    "CT": ("mixed", "Частичная переносимость"),
                    "CC": ("ancestral", "Предковый вариант - непереносимость лактозы"),
                    "AA": ("european_pastoral", "Переносит лактозу"),
                    "AG": ("mixed", "Частичная переносимость"),
                    "GG": ("ancestral", "Непереносимость лактозы"),
                }
            },
            "rs182549": {
                "gene": "LCT",
                "description": "Персистенция лактазы (арабский/африканский)",
                "interpretation": {
                    "CC": ("arabic_african", "Переносит лактозу (арабский/африканский вариант)"),
                    "CT": ("mixed", "Частичная"),
                    "TT": ("ancestral", "Предковый вариант"),
                }
            },
        }
    },

    "neanderthal": {
        "name": "Неандертальское наследие",
        "snps": {
            "rs2298813": {
                "gene": "Neanderthal",
                "description": "Направление роста волос",
                "interpretation": {
                    "AA": ("neanderthal", "Неандертальский вариант"),
                    "AG": ("mixed", "Смешанный"),
                    "GG": ("modern", "Современный человек"),
                }
            },
            "rs3917862": {
                "gene": "Neanderthal",
                "description": "Свёртываемость крови",
                "interpretation": {
                    "TT": ("neanderthal", "Неандертальский вариант"),
                    "CT": ("mixed", "Смешанный"),
                    "CC": ("modern", "Современный"),
                }
            },
        }
    },
}

# Y-хромосомные гаплогруппы (только для мужчин)
Y_HAPLOGROUPS = {
    "rs9786184": ("R1b", "M269 - Западная Европа"),
    "rs17250804": ("R1a", "M17/M198 - Восточная Европа, Азия"),
    "rs9341296": ("I1", "M253 - Скандинавия"),
    "rs2032597": ("I2", "M438 - Балканы, Восточная Европа"),
    "rs34126399": ("J1", "M267 - Ближний Восток"),
    "rs34305635": ("J2", "M172 - Средиземноморье"),
    "rs2032604": ("E1b1b", "M35 - Африка, Средиземноморье"),
    "rs2032636": ("G", "M201 - Кавказ, Европа"),
    "rs9341301": ("N", "M231 - Финно-угры, Сибирь"),
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


def analyze_ancestry(genome):
    """Analyze ancestry markers"""
    results = {}

    for category, cat_info in ANCESTRY_SNPS.items():
        cat_results = []
        for snp_id, snp_info in cat_info['snps'].items():
            result = {
                'snp_id': snp_id,
                'gene': snp_info['gene'],
                'description': snp_info['description'],
                'found': False,
                'genotype': None,
                'origin': None,
                'interpretation': None
            }

            if snp_id in genome:
                result['found'] = True
                genotype = genome[snp_id]['genotype']
                result['genotype'] = genotype

                # Try to find interpretation
                interp = snp_info.get('interpretation', {})
                for gt in [genotype, genotype[::-1] if len(genotype) == 2 else genotype]:
                    if gt in interp:
                        result['origin'], result['interpretation'] = interp[gt]
                        break

            cat_results.append(result)
        results[category] = cat_results

    return results


def analyze_y_haplogroup(genome):
    """Analyze Y-chromosome haplogroup markers"""
    results = []

    for snp_id, (haplogroup, description) in Y_HAPLOGROUPS.items():
        if snp_id in genome:
            genotype = genome[snp_id]['genotype']
            # Check if derived allele present (simplified)
            results.append({
                'snp_id': snp_id,
                'haplogroup': haplogroup,
                'description': description,
                'genotype': genotype,
                'found': True
            })

    return results


def determine_sex(genome):
    """Determine genetic sex from Y-chromosome markers"""
    # Check for Y chromosome SNPs
    y_snps = [snp for snp, data in genome.items() if data['chromosome'] == 'Y']
    return 'male' if len(y_snps) > 10 else 'female'


def generate_report(results, y_results, sex, genome):
    """Generate ancestry report"""
    report = []
    report.append("# 🌍 Анализ происхождения")
    report.append(f"\nДата анализа: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append(f"\nГенетический пол: **{'Мужской' if sex == 'male' else 'Женский'}**")
    report.append("\n---\n")

    report.append("## ⚠️ Важно\n")
    report.append("- Это упрощённый анализ на основе отдельных SNP")
    report.append("- Для точного определения происхождения нужен анализ тысяч маркеров")
    report.append("- Рекомендуется загрузить данные в GEDmatch или DNA.Land для детального анализа\n")
    report.append("---\n")

    # Ancestry Informative Markers
    report.append("## 🧬 Информативные маркеры происхождения (AIMs)\n")

    for category, cat_results in results.items():
        cat_name = ANCESTRY_SNPS[category]['name']
        report.append(f"### {cat_name}\n")
        report.append("| SNP | Ген | Генотип | Интерпретация |")
        report.append("|-----|-----|---------|---------------|")

        for r in cat_results:
            if r['found']:
                interp = r['interpretation'] or 'Нет данных'
                report.append(f"| {r['snp_id']} | {r['gene']} | **{r['genotype']}** | {interp} |")
            else:
                report.append(f"| {r['snp_id']} | {r['gene']} | - | Не найден |")
        report.append("")

    # Y-haplogroup (if male)
    if sex == 'male' and y_results:
        report.append("## 👨 Y-хромосомная гаплогруппа (отцовская линия)\n")
        report.append("| SNP | Гаплогруппа | Генотип | Описание |")
        report.append("|-----|-------------|---------|----------|")

        for r in y_results:
            report.append(f"| {r['snp_id']} | {r['haplogroup']} | {r['genotype']} | {r['description']} |")

        report.append("\n**Примечание**: Для точного определения Y-гаплогруппы рекомендуется YFull или FTDNA")
        report.append("")

    # Summary
    report.append("\n---\n")
    report.append("## 📊 Сводка\n")

    # Count origins
    origin_counts = {}
    for cat_results in results.values():
        for r in cat_results:
            if r['origin']:
                origin_counts[r['origin']] = origin_counts.get(r['origin'], 0) + 1

    origin_labels = {
        'european': '🇪🇺 Европейское',
        'northern_european': '🇸🇪 Северно-европейское',
        'european_pastoral': '🐄 Европейское пастушеское',
        'east_asian': '🇯🇵 Восточноазиатское',
        'african': '🌍 Африканское',
        'arabic_african': '🏜️ Арабское/Африканское',
        'neanderthal': '🦴 Неандертальское',
        'ancestral': '🌿 Предковое',
        'mixed': '🔀 Смешанное',
        'modern': '👤 Современный человек',
        'other': '• Другое'
    }

    if origin_counts:
        for origin, count in sorted(origin_counts.items(), key=lambda x: -x[1]):
            label = origin_labels.get(origin, origin)
            report.append(f"- {label}: {count} маркеров")

    report.append("\n---\n")
    report.append("## 🔗 Рекомендуемые ресурсы\n")
    report.append("- **GEDmatch**: gedmatch.com - бесплатное сравнение с базами данных")
    report.append("- **DNA.Land**: dna.land - бесплатный анализ происхождения")
    report.append("- **YFull**: yfull.com - детальный Y-хромосомный анализ")
    report.append("- **FTDNA**: familytreedna.com - трансфер для гаплогрупп")

    return '\n'.join(report)


def main():
    print("=" * 60)
    print("АНАЛИЗ ПРОИСХОЖДЕНИЯ ПО ГЕНОМУ 23andMe")
    print("=" * 60)

    print("\n[1/4] Загрузка генома...")
    genome = load_genome()
    print(f"      Загружено {len(genome)} SNP")

    print("\n[2/4] Определение пола...")
    sex = determine_sex(genome)
    print(f"      Генетический пол: {'мужской' if sex == 'male' else 'женский'}")

    print("\n[3/4] Анализ маркеров происхождения...")
    results = analyze_ancestry(genome)

    y_results = []
    if sex == 'male':
        print("      Анализ Y-хромосомы...")
        y_results = analyze_y_haplogroup(genome)

    print("\n[4/4] Генерация отчёта...")
    report = generate_report(results, y_results, sex, genome)

    report_path = f"{REPORTS_PATH}/ancestry/report.md"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"      → {report_path}")

    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)


if __name__ == "__main__":
    main()
