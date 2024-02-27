import re
### pandas telepítés ellenőrzése

print("pandas modul ellenőrzése...")
try:
    import pandas as pd
    print("OK!")
except ImportError:
    print("Pandas nincs telepítve. Telepítés...")
    try:
        # Pandas telepítése
        import subprocess

        subprocess.check_call(["pip", "install", "pandas"])
        print("Pandas sikeresen telepítve.")

        # Újraimportálás a frissen telepített Pandas használatához
        import pandas as pd
    except Exception as e:
        print(f"Hiba a Pandas telepítése során: {e}")
        exit()

# A sam file elérési útja
file_path = "/mnt/e/armin/MPOX_KSHV_cage_A_20230131_GCGG.Aligned.out.sam"
# A kimeneti sam file elérési útja .sam a vége!!!
out_path = "/mnt/"

def process_file(file_path):
    sam_data = pd.DataFrame()
    try:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')

                if '@' in columns[0]:
                    header = columns
                    sam_data = pd.concat([sam_data, pd.DataFrame([header])], ignore_index=True)

                else:
                    CIGAR = columns[5] if len(columns) > 5 else ''
                    SEQENCE = columns[9] if len(columns) > 5 else ''
                    QUALITY = columns[10] if len(columns) > 5 else ''
                    num_of_s = CIGAR.count('S')

                    # amikor mindkét vége softclip
                    if num_of_s == 2:
                        start_s_index = CIGAR.find('S')
                        end_s_index = CIGAR.rfind('S')
                        # Két S van, mindkét S előtti számot írja ki
                        rcutvalue = int(re.split('[A-Z]', CIGAR[start_s_index + 1:end_s_index])[1])
                        lcutvalue = int(CIGAR[:start_s_index].strip() if start_s_index > 0 else 0)
                        # A cutvalue értékekkel vágni kell
                        corrected_seq = SEQENCE[lcutvalue:len(SEQENCE) - rcutvalue]
                        corrected_quality = QUALITY[lcutvalue:len(QUALITY) - rcutvalue]
                        columns[5] = CIGAR[len(str(lcutvalue))+1:len(CIGAR) - (len(str(lcutvalue)) + 1)]
                        columns[9] = corrected_seq
                        columns[10] = corrected_quality
                        sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)
                    #amikor csak az egyik vége softclip
                    elif num_of_s == 1:
                        start_s_index = CIGAR.find('S')
                        start_m_index = CIGAR.rfind('M')

                        # amikor a softclip a jobb végén van
                        if start_s_index > start_m_index:
                            rcutvalue = int(re.split('[M]', CIGAR[start_m_index + 1:start_s_index])[0])
                            lcutvalue = 0
                            corrected_seq = SEQENCE[lcutvalue:len(SEQENCE) - rcutvalue]
                            corrected_quality = QUALITY[lcutvalue:len(QUALITY) - rcutvalue]
                            columns[5] = CIGAR[0:len(CIGAR) - (len(str(lcutvalue)) + 1)]
                            columns[9] = corrected_seq
                            columns[10] = corrected_quality
                            sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)

                        # amikor a softclip a bal végén van
                        else:
                            rcutvalue = 0
                            lcutvalue = int(CIGAR[:start_s_index].strip())
                            corrected_seq = SEQENCE[lcutvalue:len(SEQENCE) - rcutvalue]
                            corrected_quality = QUALITY[lcutvalue:len(QUALITY) - rcutvalue]
                            columns[5] = CIGAR[len(str(lcutvalue))+1:len(CIGAR)]
                            columns[9] = corrected_seq
                            columns[10] = corrected_quality
                            sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)


                    # amikor a read unmapped ezért nem változtatunk az adaton
                    elif num_of_s == 0:
                        sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)
    except FileNotFoundError:
        print(f"A fájl nem található: {file_path}")
    except Exception as e:
        print(f"Hiba történt: {e}")
    print(sam_data)

# Hívja meg a függvényt a fájl feldolgozásához
output = process_file(file_path)

output.to_csv(out_path, sep='\t', index=False, header=False)