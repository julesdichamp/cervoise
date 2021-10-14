import matplotlib.pyplot as plt

betaTab = [1, 5, 10]
varTab1 = [469.0, 590.0, 657.0]
varTab2 = [475.0, 578.0, 650.0]
varTab3 = [484.0, 562.0, 646.0]

plt.figure(figsize=(16, 10))
plt.plot(betaTab, varTab1, label = "beta/beta* = 1")
plt.plot(betaTab, varTab2, label = "beta/beta* = 2")
plt.plot(betaTab, varTab3, label = "beta/beta* = 3")
plt.title("Parametric study, S_ini/Sss_ini = 10")
plt.legend()
plt.xlabel("Ratio k*/k")
plt.ylabel("Max nb of infected people")
plt.tight_layout()
plt.show()
