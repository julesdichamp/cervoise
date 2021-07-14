import matplotlib.pyplot as plt

betaTab = [1, 2, 3]
varTab = [29.40260768269755, 16.25438868927738, 8.147890516279395]

plt.plot(betaTab, varTab)
plt.xlabel("Values of Beta")
plt.ylabel("avg return time")
plt.tight_layout()
plt.show()
