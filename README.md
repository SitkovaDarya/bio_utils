# bio_utils

Небольшой набор биоинформатических утилит для работы с файлами форматов
FASTA, FASTQ, BLAST и GenBank.

## Установка

```bash
git clone https://github.com/SitkovaDarya/bio_utils.git
cd bio_utils
```

## Использование

```python
from bio_utils import (
    filter_fastq,
    parse_blast_output,
    select_genes_from_gbk_to_fasta,
)

# Фильтрация FASTQ по GC-составу, длине и среднему качеству
filter_fastq(
    "input.fastq",
    "out.fastq",
    gc_bounds=(20, 80),
    length_bounds=(50, 200),
    quality_threshold=20,
)
# Результат сохранится в filtered/out.fastq

# Извлечь лучший хит BLAST по каждому запросу
parse_blast_output("blast_results.txt", "hits.txt")

# Извлечь соседние гены вокруг заданного гена из GBK-аннотации
select_genes_from_gbk_to_fasta(
    "genome.gbk", "thrC", n_before=1, n_after=2, output_fasta="neighbors.fasta"
)
```

## Возможности

- **`filter_fastq`** — фильтрация ридов FASTQ по трём критериям:
  - GC-состав (в процентах);
  - длина;
  - среднее качество (шкала Phred33).

  Работает «на лету»: в памяти держится только один рид за раз.

- **`parse_blast_output`** — читает текстовый отчёт BLAST, для каждого
  запроса берёт название лучшего совпадения (первая строка столбца
  `Description`) и сохраняет все результаты одним столбцом,
  отсортированным по алфавиту.

- **`select_genes_from_gbk_to_fasta`** — читает аннотацию генома в формате
  GenBank, находит гены интереса и сохраняет в FASTA белковые
  последовательности соседних генов (N генов до и M после).

## Структура пакета

- `bio_utils.py` — точка входа, публичные функции
- `helpers/fastq.py` — вспомогательные функции для FASTQ-фильтрации
- `helpers/bio_files.py` — вспомогательные функции для парсинга FASTA/BLAST/GBK

## Требования

- Python 3.10+
- Только стандартная библиотека (внешние зависимости не нужны)

## Контакты

SitkovaDarya — https://github.com/SitkovaDarya