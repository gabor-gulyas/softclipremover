import re
### pandas installation verification

print("pandas module check...")
try:
    import pandas as pd
    print("OK!")
except ImportError:
    print("Pandas are not installed. Installation...")
    try:
        # Pandas telepítése
        import subprocess

        subprocess.check_call(["pip", "install", "pandas"])
        print("Pandas successfully installed.")

        # Re-import to use freshly installed Pandas
        import pandas as pd
    except Exception as e:
        print(f"Error during the installation of Pandas: {e}")
        exit()

# input sam file path
file_path = "/mnt/e/armin/MPOX_KSHV_cage_A_20230131_GCGG.Aligned.out.sam"
# output path. if u would like sam file, dont forget the sam extension
out_path = "/mnt/"

def softclipremover(file_path):
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

                    # Softclip are in both end
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

                        # Soft clip is right end
                        if start_s_index > start_m_index:
                            rcutvalue = int(re.split('[M]', CIGAR[start_m_index + 1:start_s_index])[0])
                            lcutvalue = 0
                            corrected_seq = SEQENCE[lcutvalue:len(SEQENCE) - rcutvalue]
                            corrected_quality = QUALITY[lcutvalue:len(QUALITY) - rcutvalue]
                            columns[5] = CIGAR[0:len(CIGAR) - (len(str(lcutvalue)) + 1)]
                            columns[9] = corrected_seq
                            columns[10] = corrected_quality
                            sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)

                        # Softclip in left end
                        else:
                            rcutvalue = 0
                            lcutvalue = int(CIGAR[:start_s_index].strip())
                            corrected_seq = SEQENCE[lcutvalue:len(SEQENCE) - rcutvalue]
                            corrected_quality = QUALITY[lcutvalue:len(QUALITY) - rcutvalue]
                            columns[5] = CIGAR[len(str(lcutvalue))+1:len(CIGAR)]
                            columns[9] = corrected_seq
                            columns[10] = corrected_quality
                            sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)


                    # When the reads are unmapped we dont change the data
                    elif num_of_s == 0:
                        sam_data = pd.concat([sam_data, pd.DataFrame([columns])], ignore_index=True)
    except FileNotFoundError:
        print(f"The file is not found {file_path}")
    except Exception as e:
        print(f"Error: {e}")
    return sam_data

# process
output = softclipremover(file_path)

output.to_csv(out_path, sep='\t', index=False, header=False)
