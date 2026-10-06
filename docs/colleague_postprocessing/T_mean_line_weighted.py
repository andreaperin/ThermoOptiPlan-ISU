import pandas as pd
import os

# Ordnerpfad anpassen, in dem sich die Excel-Dateien befinden
input_folder = "C:/Users/niehage/Documents/01_Projekt/01_Modell/03_MehrschichtModell_line" 

# Dateinamen (du kannst beliebig viele ergänzen)
files = ["Tm_Multi_BW_line_200_S2.csv", "Tm_Multi_BW_line_200_S3.csv", "Tm_Multi_BW_line_200_S4.csv"]

# Gewichte entsprechend der Fließraten
weights = [0.245, 0.43, 0.325]

# Alle Dateien einlesen
dfs = []
for i, fname in enumerate(files):
    file_path = os.path.join(input_folder, fname)
    df = pd.read_csv(file_path)
    df = df.rename(columns={'T_mean': f'T_mean_{i+1}'})
    dfs.append(df)

# Merge über 'Time'
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on='Time')

# Gewichtete Mittelung
t_mean_cols = [f'T_mean_{i+1}' for i in range(len(weights))]
merged['T_mean_avg'] = sum(
    w * merged[col] for w, col in zip(weights, t_mean_cols)
)

# Ausgabe
output_path = os.path.join(
    input_folder, "T_mean_Multi_BW_line_200_weighted.csv"
)
merged[['Time', 'T_mean_avg']].to_csv(output_path, index=False)

print(f"Gewichtete Mittelwert-Datei gespeichert unter: {output_path}")
