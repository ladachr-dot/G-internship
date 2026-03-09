import pandas as pd
import numpy as np
import os
import re
import requests
import time
import math
from tabulate import tabulate
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
import csv
import pandas as pd
import gspread
from google.oauth2.service_account import Credentials
# ------------------ Helpers ------------------

def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns: return c
    return None


_float_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def parse_first_float(x):
    if x is None or (isinstance(x, float) and np.isnan(x)): return np.nan
    s = str(x).strip()
    if s.upper() in ("", "-", "NA", "NAN", "NONE", "NR"): return np.nan
    m = _float_re.search(s)
    return float(m.group(0)) if m else np.nan


def parse_beta_signed(x):
    if x is None or (isinstance(x, float) and np.isnan(x)): return np.nan
    s = str(x).strip().lower()
    m = _float_re.search(s)
    if not m: return np.nan
    v = float(m.group(0))
    if "decrease" in s and v > 0: v = -v
    return v


def extract_rsid(x):
    if x is None: return None
    m = re.search(r"(rs\d+)", str(x), flags=re.IGNORECASE)
    return m.group(1).lower() if m else None


# ------------------ API & Processing ------------------

def fetch_variant_data(rsid):
    server = "https://rest.ensembl.org"
    ext = f"/variation/human/{rsid}?pops=1"
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if r.status_code == 429:
            time.sleep(1)
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        return r.json() if r.ok else None
    except:
        return None


def build_enriched_table(file_path, top_n=20):
    file_path = file_path.strip().replace('"', '').replace("'", "")
    search_paths = [file_path, os.path.join(os.path.expanduser("~"), "Downloads", file_path)]
    found_path = next((p for p in search_paths if os.path.exists(p)), None)

    if not found_path:
        raise FileNotFoundError(f"Файл '{file_path}' не найден.")

    df = pd.read_csv(found_path, sep="\t", low_memory=False)

    var_col = pick_col(df, ["riskAllele", "SNPS", "variantId"])
    or_col = pick_col(df, ["orValue", "OR", "odds_ratio"])
    beta_col = pick_col(df, ["beta", "BETA", "Beta"])
    raf_col = pick_col(df, ["riskFrequency", "RAF", "raf"])

    df["raf_num"] = df[raf_col].apply(parse_first_float)
    df = df[df["raf_num"].notna()].copy()

    df["or_num"] = df[or_col].apply(parse_first_float) if or_col else np.nan
    df["beta_num"] = df[beta_col].apply(parse_beta_signed) if beta_col else np.nan

    or_count = df["or_num"].notna().sum()
    beta_count = df["beta_num"].notna().sum()
    effect_used = "OR" if or_count >= beta_count else "BETA"

    if effect_used == "OR":
        df["score"] = df["or_num"].apply(lambda v: abs(math.log(v)) if (v and v > 0) else np.nan)
    else:
        df["score"] = df["beta_num"].abs()

    df["rsid"] = df[var_col].apply(extract_rsid)

    df_filtered = (
        df
        .dropna(subset=["rsid"])
        .drop_duplicates(subset=["rsid"])
    )

    return effect_used, df_filtered.head(top_n).copy()


# ------------------ Main Pipeline ------------------

if __name__ == "__main__":
    file_input = input("Введите имя TSV файла: ").strip()
    top_n_in = input("Сколько вариантов в таблицу (по умолчанию 20): ").strip()
    top_n = int(top_n_in) if top_n_in else 20

    try:
        effect_used, final_df = build_enriched_table(file_input, top_n=top_n)
        print(f"\n[!] Используем эффект {effect_used}. Начинаю сбор данных из Ensembl...")

        api_rows = []
        euro_codes = {'CEU', 'FIN', 'GBR', 'IBS', 'TSI'}

        for _, row in final_df.iterrows():
            rsid = row["rsid"]
            print(f" -> Обработка {rsid}...")

            data = fetch_variant_data(rsid)
            # Заготовка строки
            entry = {
                "consequence": "Error/Not Found",
                "clin_sig": "N/A"
            }
            # Заготовка под частоты (по умолчанию N/A)
            for code in euro_codes: entry[f"freq_{code}"] = "N/A"

            if data:
                entry["consequence"] = data.get('most_severe_consequence', 'N/A')
                entry["clin_sig"] = ", ".join(data.get('clinical_significance', [])) or 'N/A'

                # Вытаскиваем частоты
                for pop in data.get('populations', []):
                    pop_name = pop.get('population', '')
                    code = pop_name.split(':')[-1] if ':' in pop_name else pop_name
                    if code in euro_codes:
                        # Записываем частоту и аллель для наглядности
                        entry[f"freq_{code}"] = f"{pop.get('frequency')} ({pop.get('allele')})"

            api_rows.append(entry)

        # Склеиваем всё в одну большую таблицу
        api_df = pd.DataFrame(api_rows)
        final_df = pd.concat([final_df.reset_index(drop=True), api_df], axis=1)

        # Формируем список колонок для вывода на экран (укороченный)
        display_cols = ["rsid", "raf_num", "score", "consequence", "freq_CEU", "freq_GBR"]
        print("\n" + "=" * 80)
        print("ТАБЛИЦА")
        print("=" * 80)
        print(tabulate(final_df[display_cols].head(10), headers='keys', tablefmt='psql', showindex=False))

        # Сохраняем ОДИН файл
        output_name = "genomics_report.csv"
        final_df.to_csv(output_name, index=False)

        # Сохраняем запрос для NCBI
        with open("ncbi_query.txt", "w") as f:
            f.write(" OR ".join(final_df["rsid"].tolist()))

        print(f"\n[+] ГОТОВО! Всё сохранено в один файл: {output_name}")
        print(f"[+] NCBI query сохранен в: ncbi_query.txt")



 # Перенос в гугл таблицы
        CSV_FILE = 'genomics_report.csv'  # Путь к CSV файлу
        CREDENTIALS_FILE = 'credentials.json'  # Файл с ключами от сервисного аккаунта
        SPREADSHEET_NAME = 'Gen'  # Имя таблицы на Google Диске
        SHEET_NAME = 'Лист1'  # Имя листа внутри таблицы

        # 1. Авторизация и подключение к Google Sheets API
        scope = ['https://www.googleapis.com/auth/spreadsheets',
                 'https://www.googleapis.com/auth/drive']
        creds = Credentials.from_service_account_file(CREDENTIALS_FILE, scopes=scope)
        client = gspread.authorize(creds)

        # 2. Открываем таблицу по имени
        try:
            spreadsheet = client.open(SPREADSHEET_NAME)
        except gspread.SpreadsheetNotFound:
            print(
                f"Ошибка: Таблица с именем '{SPREADSHEET_NAME}' не найдена. Убедитесь, что вы поделились ею с сервисным аккаунтом.")
            exit()

        # 3. Выбираем лист
        try:
            worksheet = spreadsheet.worksheet(SHEET_NAME)
        except gspread.WorksheetNotFound:
            print(f"Лист '{SHEET_NAME}' не найден, создаем новый.")
            worksheet = spreadsheet.add_worksheet(title=SHEET_NAME, rows=100, cols=20)

        # 4. Читаем CSV файл с помощью pandas
        try:
            df = pd.read_csv(CSV_FILE, encoding='utf-8')  # При необходимости измените кодировку

            # Заменяем NaN на пустые строки (иначе он посылает далеко)
            df = df.fillna('')

            print(f"CSV файл '{CSV_FILE}' прочитан. Найдено {len(df)} строк.")
            print(f"Все NaN значения заменены на пустые строки.")

        except FileNotFoundError:
            print(f"Ошибка: CSV файл '{CSV_FILE}' не найден.")
            exit()
        except Exception as e:
            print(f"Ошибка при чтении CSV файла: {e}")
            exit()

        # 5. Преобразуем DataFrame в формат для загрузки (список списков)
        #    Первый список - заголовки, остальные - данные
        data_to_upload = [df.columns.values.tolist()] + df.values.tolist()

        # 6. Очищаем лист и загружаем новые данные
        try:
            worksheet.clear()  # Очищаем лист от старых данных
            worksheet.update(data_to_upload, 'A1')  # Сначала значения, потом диапазон
            print(
                f"Данные из '{CSV_FILE}' успешно загружены в Google Таблицу '{SPREADSHEET_NAME}', лист '{SHEET_NAME}'.")
        except Exception as e:
            print(f"Ошибка при загрузке данных в Google Таблицу: {e}")
            exit()
    except Exception as e:
        print(f"\n[!] ОШИБКА: {e}")
