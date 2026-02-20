import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def apply_grey_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='0.92')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("0.85")

def main():
	data = pd.read_excel("data_processed.xlsx")
	time = data["Time (hours)"]
	epsilon_1 = data["Epsilon_1 (%)"]
	epsilon_3 = data["Epsilon_3 (%)"]
	sigma_1 = data["Sigma_1 (MPa)"]
	sigma_3 = data["Sigma_3 (MPa)"]

	fig, ax = plt.subplots(1, 1, figsize=(8, 4))
	fig.subplots_adjust(top=0.880, bottom=0.115, left=0.075, right=0.900, hspace=0.245, wspace=0.235)

	ax.plot(time[-2:-1], epsilon_1[-2:-1], '-', color="steelblue", linewidth=2.0, label="Axial strain")
	ax.plot(time[-2:-1], epsilon_3[-2:-1], '-', color="lightcoral", linewidth=2.0, label="Radial strain")
	ax.plot(time, sigma_1, '-', color="0.5", linewidth=2.0, label="Axial stress")
	ax.plot(time, sigma_3, '-', color="0.2", linewidth=2.0, label="Radial stress")
	ax.set_xlabel("Time (hours)", size=12, fontname="serif")
	ax.set_ylabel("Stress (MPa)", size=12, fontname="serif")
	ax.legend(loc=0, shadow=True, fancybox=True, bbox_to_anchor=[0.95, 1.17], ncol=4)
	ax.grid(False)

	axt = ax.twinx()

	axt.plot(time, epsilon_1, '-', color="steelblue", linewidth=2.0)
	axt.plot(time, epsilon_3, '-', color="lightcoral", linewidth=2.0)
	axt.set_ylabel("Total strain (%)", size=12, fontname="serif")

	apply_grey_theme(fig, [ax], transparent=True)

	plt.show()

if __name__ == "__main__":
	main()