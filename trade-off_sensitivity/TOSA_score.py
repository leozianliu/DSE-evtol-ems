import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Reliability scores
rs_concept_1 = np.array([5, 4, 4, 4, 5, 5, 4])
rs_concept_2 = np.array([4, 3, 3, 3, 3, 3, 3])
rs_concept_3 = np.array([3, 3, 1, 4, 2, 4, 3])

#Maintanability scores
ms_concept_1 = np.array([4, 4, 3, 5, 4, 3, 3])
ms_concept_2 = np.array([4, 2, 3, 3, 2, 3, 4])
ms_concept_3 = np.array([2, 2, 2, 2, 2, 4, 3])

#average and std scores for realibility
rs_average_1 = np.mean(rs_concept_1.astype(float))
rs_std_1 = np.std(rs_concept_1.astype(float))
rs_average_2 = np.mean(rs_concept_2.astype(float))
rs_std_2 = np.std(rs_concept_2.astype(float))
rs_average_3 = np.mean(rs_concept_3.astype(float))
rs_std_3 = np.std(rs_concept_3.astype(float))

#average and std scores for maintanability
ms_average_1 = np.mean(ms_concept_1.astype(float))
ms_std_1 = np.std(ms_concept_1.astype(float))
ms_average_2 = np.mean(ms_concept_2.astype(float))
ms_std_2 = np.std(ms_concept_2.astype(float))
ms_average_3 = np.mean(ms_concept_3.astype(float))
ms_std_3 = np.std(ms_concept_3.astype(float))

print(rs_average_1, rs_std_1, rs_average_2, rs_std_2, rs_average_3, rs_std_3)
print(ms_average_1, ms_std_1, ms_average_2, ms_std_2, ms_average_3, ms_std_3)

criteria_weights = np.array([3., 1., 2., 2., 2.])
#reliability, maintainability, controllability, peak power, takeoff weight
w_r = criteria_weights[0]
w_m = criteria_weights[1]

#propagation of uncertainty (full concept)
wr = (w_r)/(sum(criteria_weights))
wm = (w_m)/(sum(criteria_weights))

unc_1 = np.sqrt((wr**2)*(rs_std_1**2) + (wm**2)*(ms_std_1**2))
unc_2 = np.sqrt((wr**2)*(rs_std_2**2) + (wm**2)*(ms_std_2**2))
unc_3 = np.sqrt((wr**2)*(rs_std_3**2) + (wm**2)*(ms_std_3**2))

#final scores from trade-off (no rounding)
concept_1 = np.array([float(rs_average_1), float(ms_average_1), 5, 1, 2])
concept_2 = np.array([rs_average_2, ms_average_2, 4, 3, 4])
concept_3 = np.array([rs_average_3, ms_average_3, 3, 5, 4])


#averaged scores (these values will vary from the ones in the trade-off, since they are not rounded)
concept_1_final = float(sum(concept_1 * criteria_weights)) / float(sum(criteria_weights))
concept_2_final = float(sum(concept_2 * criteria_weights)) / float(sum(criteria_weights))
concept_3_final = float(sum(concept_3 * criteria_weights)) / float(sum(criteria_weights))

#final scores with uncertainty
print("Concept 1", concept_1_final, "+/-", unc_1)
print("Concept 2", concept_2_final, "+/-", unc_2)
print("Concept 3", concept_3_final, "+/-", unc_3)


data_1 = np.array([concept_1_final + unc_1, concept_1_final, concept_1_final - unc_1])
data_2 = np.array([concept_2_final + unc_2, concept_2_final, concept_2_final - unc_2])
data_3 = np.array([concept_3_final + unc_3, concept_3_final, concept_3_final - unc_3])

x = [1, 2, 3]
y = [concept_1_final, concept_2_final, concept_3_final]
yerr = [unc_1, unc_2, unc_3]
y1 = concept_3_final - unc_3
y2 = concept_3_final 

plt.figure(figsize=(8, 5))
plt.errorbar(x, y, yerr=yerr, fmt='x', capsize=5, capthick=2, elinewidth=2)
plt.axhline(y=y1, linestyle='dashed', color='lightblue')
plt.axhline(y=y2, linestyle='dashed', color='black')
plt.xticks(x, ['Concept 1', 'Concept 2', 'Concept 3'])
plt.ylabel('Final Score')
plt.show()