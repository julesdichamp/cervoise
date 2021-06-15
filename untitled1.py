import matplotlib.pyplot as plt

timeTab = [180, 1800]
runningTime1 = [31.28297766049703,261.39304971694946]
runningTime2 = [32.58021652698517,261.7156163454056]
runningTime3 = [31.509316086769104,269.0880151987076]


plt.plot(timeTab, runningTime1, "ro",label="full changed")
plt.plot(timeTab, runningTime2, "go", label="only dichotomy")
plt.plot(timeTab, runningTime2, "bo",label="no change")
plt.title("config_VALIDATION_SIR_compartments_fullConnectedRegions_agent_based.ini, t_end=50")
plt.xlabel("S_ini")
plt.legend()

plt.show()
