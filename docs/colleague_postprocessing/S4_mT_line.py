from paraview.simple import *
import numpy as np
import csv
import os

# === Parameter ===
mesh_path = "C:/Users/niehage/Documents/01_Projekt/01_Modell/03_MehrschichtModell_line/Multi_BW_line_200.pvd"
output_csv = "C:/Users/niehage/Documents/01_Projekt/01_Modell/03_MehrschichtModell_line/Tm_Multi_BW_line_200_S4.csv"
fieldname = 'T'

# Linienendpunkte
p1 = np.array([2250.0, 0.0, -1364.6])
p2 = np.array([2250.0, 0.0, -1374.8])

# Toleranz für "auf der Linie"
eps = 1e-4   # ggf. an Meshgröße anpassen

z_top = -1364.6
z_bottom = -1374.8

# =========================
# Output-Pfad
# =========================

os.makedirs(os.path.dirname(output_csv), exist_ok=True)

# =========================
# Mesh & Zeit
# =========================

mesh = OpenDataFile(mesh_path)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

timekeeper = GetTimeKeeper()
times = timekeeper.TimestepValues

rows = []

# Richtungsvektor der Linie
line_vec = p2 - p1
line_len2 = np.dot(line_vec, line_vec)

# =========================
# Zeitschleife
# =========================

for t in times:
    scene.AnimationTime = t
    mesh.UpdatePipeline(time=t)

    data = servermanager.Fetch(mesh)
    points = data.GetPoints()
    pdata = data.GetPointData()

    if not pdata.HasArray(fieldname):
        rows.append({'Time': t, 'T_mean': None})
        continue

    arr = pdata.GetArray(fieldname)

    values = []
    n_total = points.GetNumberOfPoints()

    for i in range(n_total):
        P = np.array(points.GetPoint(i))
        z = P[2]

        if not (z_bottom <= z <= z_top):
            continue

        # Projektion auf Linie
        t_proj = np.dot(P - p1, line_vec) / line_len2
        t_proj = max(0.0, min(1.0, t_proj))

        closest = p1 + t_proj * line_vec
        dist = np.linalg.norm(P - closest)

        if dist <= eps:
            values.append(arr.GetValue(i))

    if len(values) == 0:
        print(f"⚠️ t={t:.2f}: 0 Linien-Knoten")
        rows.append({'Time': t, 'T_mean': None})
    else:
        mean_val = np.mean(values)
        print(f"t={t:.2f} → {len(values)} Knoten → ⌀T = {mean_val:.3f}")
        rows.append({'Time': t, 'T_mean': mean_val})

# =========================
# CSV schreiben
# =========================

with open(output_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['Time', 'T_mean'])
    writer.writeheader()
    writer.writerows(rows)

print("✅ CSV gespeichert:", output_csv)
