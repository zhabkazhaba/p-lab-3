import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../data/out_times.csv")

plt.figure(figsize=(10, 6))

plt.plot(df["Size"], df["CubeCG"], marker='o', label='CubeCG')
plt.plot(df["Size"], df["LCGShuffle"], marker='s', label='LCGShuffle')
plt.plot(df["Size"], df["QuickBBS"], marker='^', label='QuickBBS')
plt.plot(df["Size"], df["MT19937"], marker='x', linestyle='--', label='std::mt19937')

plt.title("Сравнение скорости генерации", fontsize=14)
plt.xlabel("Объем выборки", fontsize=12)
plt.ylabel("Время генерации (мс)", fontsize=12)

plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=10)

plt.tight_layout()

plt.savefig("../data/rando_plot.png", dpi=300)
plt.show()