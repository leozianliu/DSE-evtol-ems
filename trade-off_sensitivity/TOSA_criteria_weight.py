import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

#weights
to_weights = np.array([3, 1, 2, 2, 2])
#changing the first weight
to_weights_2 = np.array([2, 1, 2, 2, 2])
to_weights_3 = np.array([4, 1, 2, 2, 2])
#changing the second weight
to_weights_4 = np.array([3, 2, 2, 2, 2])
#changing the third weight
to_weights_5 = np.array([3, 1, 1, 2, 2])
to_weights_6 = np.array([3, 1, 3, 2, 2])
#changing the fourth weight
to_weights_7 = np.array([3, 1, 2, 1, 2])
to_weights_8 = np.array([3, 1, 2, 3, 2])
#changing the fifth weight
to_weights_9 = np.array([3, 1, 2, 2, 1])
to_weights_10 = np.array([3, 1, 2, 2, 3])


#total weights
total_weight = np.sum(to_weights)
total_weight_2 = np.sum(to_weights_2)
total_weight_3 = np.sum(to_weights_3)
total_weight_4 = np.sum(to_weights_4)
total_weight_5 = np.sum(to_weights_5)
total_weight_6 = np.sum(to_weights_6)
total_weight_7 = np.sum(to_weights_7)
total_weight_8 = np.sum(to_weights_8)
total_weight_9 = np.sum(to_weights_9)
total_weight_10 = np.sum(to_weights_10)

#in percentage
percentage_weights = to_weights / total_weight
percentage_weights_2 = to_weights_2 / total_weight_2
percentage_weights_3 = to_weights_3 / total_weight_3  
percentage_weights_4 = to_weights_4 / total_weight_4 
percentage_weights_5 = to_weights_5 / total_weight_5
percentage_weights_6 = to_weights_6 / total_weight_6
percentage_weights_7 = to_weights_7 / total_weight_7
percentage_weights_8 = to_weights_8 / total_weight_8
percentage_weights_9 = to_weights_9 / total_weight_9
percentage_weights_10 = to_weights_10 / total_weight_10


#scores from trade-off (do not change these)
concept_1_scores = np.array([4, 4, 5, 1, 2])
concept_2_scores = np.array([3, 3, 4, 3, 4])
concept_3_scores = np.array([3, 2, 3, 5, 4])
concept_4_scores = np.array([3, 3, 5, 2, 3])
concept_5_scores = np.array([2, 2, 3, 2, 3])

#first weights
concept_1_final = np.sum(percentage_weights * concept_1_scores)
concept_2_final = np.sum(percentage_weights * concept_2_scores)
concept_3_final = np.sum(percentage_weights * concept_3_scores)
concept_4_final = np.sum(percentage_weights * concept_4_scores)
concept_5_final = np.sum(percentage_weights * concept_5_scores)

#changing the first weight
#second weights
concept_1_final_2 = np.sum(percentage_weights_2 * concept_1_scores)
concept_2_final_2 = np.sum(percentage_weights_2 * concept_2_scores)
concept_3_final_2 = np.sum(percentage_weights_2 * concept_3_scores)
concept_4_final_2 = np.sum(percentage_weights_2 * concept_4_scores)
concept_5_final_2 = np.sum(percentage_weights_2 * concept_5_scores)

#third weights
concept_1_final_3 = np.sum(percentage_weights_3 * concept_1_scores)
concept_2_final_3 = np.sum(percentage_weights_3 * concept_2_scores)
concept_3_final_3 = np.sum(percentage_weights_3 * concept_3_scores)
concept_4_final_3 = np.sum(percentage_weights_3 * concept_4_scores)
concept_5_final_3 = np.sum(percentage_weights_3 * concept_5_scores)


#changing the second weight
#fourth weights 
concept_1_final_4 = np.sum(percentage_weights_4 * concept_1_scores)
concept_2_final_4 = np.sum(percentage_weights_4 * concept_2_scores)
concept_3_final_4 = np.sum(percentage_weights_4 * concept_3_scores)
concept_4_final_4 = np.sum(percentage_weights_4 * concept_4_scores)
concept_5_final_4 = np.sum(percentage_weights_4 * concept_5_scores)


#changing the third weight
#fifth weights
concept_1_final_5 = np.sum(percentage_weights_5 * concept_1_scores)
concept_2_final_5 = np.sum(percentage_weights_5 * concept_2_scores)
concept_3_final_5 = np.sum(percentage_weights_5 * concept_3_scores)
concept_4_final_5 = np.sum(percentage_weights_5 * concept_4_scores)
concept_5_final_5 = np.sum(percentage_weights_5 * concept_5_scores)

#sixth weights
concept_1_final_6 = np.sum(percentage_weights_6 * concept_1_scores)
concept_2_final_6 = np.sum(percentage_weights_6 * concept_2_scores)
concept_3_final_6 = np.sum(percentage_weights_6 * concept_3_scores)
concept_4_final_6 = np.sum(percentage_weights_6 * concept_4_scores)
concept_5_final_6 = np.sum(percentage_weights_6 * concept_5_scores)


#changing the fourth weight 
#seventh weights
concept_1_final_7 = np.sum(percentage_weights_7 * concept_1_scores)
concept_2_final_7 = np.sum(percentage_weights_7 * concept_2_scores)
concept_3_final_7 = np.sum(percentage_weights_7 * concept_3_scores)
concept_4_final_7 = np.sum(percentage_weights_7 * concept_4_scores)
concept_5_final_7 = np.sum(percentage_weights_7 * concept_5_scores)

#eighth weights
concept_1_final_8 = np.sum(percentage_weights_8 * concept_1_scores)
concept_2_final_8 = np.sum(percentage_weights_8 * concept_2_scores)
concept_3_final_8 = np.sum(percentage_weights_8 * concept_3_scores)
concept_4_final_8 = np.sum(percentage_weights_8 * concept_4_scores)
concept_5_final_8 = np.sum(percentage_weights_8 * concept_5_scores)


#changing the fifth weight
#ninth weights
concept_1_final_9 = np.sum(percentage_weights_9 * concept_1_scores)
concept_2_final_9 = np.sum(percentage_weights_9 * concept_2_scores)
concept_3_final_9 = np.sum(percentage_weights_9 * concept_3_scores)
concept_4_final_9 = np.sum(percentage_weights_9 * concept_4_scores)
concept_5_final_9 = np.sum(percentage_weights_9 * concept_5_scores)

#tenth weights
concept_1_final_10 = np.sum(percentage_weights_10 * concept_1_scores)
concept_2_final_10 = np.sum(percentage_weights_10 * concept_2_scores)
concept_3_final_10 = np.sum(percentage_weights_10 * concept_3_scores)
concept_4_final_10 = np.sum(percentage_weights_10 * concept_4_scores)
concept_5_final_10 = np.sum(percentage_weights_10 * concept_5_scores)


### plotting the results for visualization
#Reliability
plt.figure(figsize=(6,4))
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final, concept_2_final, concept_3_final], marker='x', label='Original Weights')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_2, concept_2_final_2, concept_3_final_2], marker='x', label='Reliability -1')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_3, concept_2_final_3, concept_3_final_3], marker='x', label='Reliability +1')
plt.ylabel('Scores')
plt.title('Reliability')
plt.legend()
plt.show()
plt.close()

#Maintainability
plt.figure(figsize=(6,4))
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final, concept_2_final, concept_3_final], marker='x', label='Original Weights')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_4, concept_2_final_4, concept_3_final_4], marker='x', label='Maintainability +1')
plt.ylabel('Scores')
plt.legend()
plt.show()
plt.close()

#Controllability
plt.figure(figsize=(6,4))
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final, concept_2_final, concept_3_final], marker='x', label='Original Weights')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_5, concept_2_final_5, concept_3_final_5], marker='x', label='Controllability -1')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_6, concept_2_final_6, concept_3_final_6], marker='x', label='Controllability +1')
plt.ylabel('Scores')
plt.legend()
plt.show()
plt.close()


#Peak power
plt.figure(figsize=(6, 4))
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final, concept_2_final, concept_3_final], marker='x', label='Original Weights')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_7, concept_2_final_7, concept_3_final_7], marker='x', label='Peak power -1')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_8, concept_2_final_8, concept_3_final_8], marker='x', label='Peak power +1')
plt.ylabel('Scores')
plt.legend()
plt.show()
plt.close()


#Takeoff weight
plt.figure(figsize=(6,4))
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final, concept_2_final, concept_3_final], marker='x', label='Original Weights')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_9, concept_2_final_9, concept_3_final_9], marker='x', label='Takeoff weight -1')
plt.scatter(['Concept 1', 'Concept 2', 'Concept 3'], [concept_1_final_10, concept_2_final_10, concept_3_final_10], marker='x', label='Takeoff weight +1')
plt.ylabel('Scores')
plt.title('Takeoff Weight')
plt.legend()
plt.show()
plt.close()


#adding to csv file 
row_0 = ['Original Weights', concept_1_final, concept_2_final, concept_3_final, concept_4_final, concept_5_final]
row_1 = ['Reliability -1', concept_1_final_2, concept_2_final_2, concept_3_final_2, concept_4_final_2, concept_5_final_2]
row_2 = ['Reliability +1', concept_1_final_3, concept_2_final_3, concept_3_final_3, concept_4_final_3, concept_5_final_3]
row_3 = ['Maintainability +1', concept_1_final_4, concept_2_final_4, concept_3_final_4, concept_4_final_4, concept_5_final_4]
row_4 = ['Controllability -1', concept_1_final_5, concept_2_final_5, concept_3_final_5, concept_4_final_5, concept_5_final_5] 
row_5 = ['Controllability +1', concept_1_final_6, concept_2_final_6, concept_3_final_6, concept_4_final_6, concept_5_final_6]
row_6 = ['Peak power -1', concept_1_final_7, concept_2_final_7, concept_3_final_7, concept_4_final_7, concept_5_final_7]
row_7 = ['Peak power +1', concept_1_final_8, concept_2_final_8, concept_3_final_8, concept_4_final_8, concept_5_final_8]
row_8 = ['Takeoff weight -1', concept_1_final_9, concept_2_final_9, concept_3_final_9, concept_4_final_9, concept_5_final_9]
row_9 = ['Takeoff weight +1', concept_1_final_10, concept_2_final_10, concept_3_final_10, concept_4_final_10, concept_5_final_10]

rows = [row_0, row_1, row_2, row_3, row_4, row_5, row_6, row_7, row_8, row_9]

#making a csv file with all the data

'''
with open('trade-off_sensitivity/tosa.csv', 'a', newline='') as file:
    writer = csv.writer(file)
    for i in range(len(rows)):
        row = rows[i]
        writer.writerow(row)
'''

#reading the csv file
df = pd.read_csv('trade-off_sensitivity/tosa.csv')
concept_1 = df.iloc[1:, 1]
concept_2 = df.iloc[1:, 2]
concept_3 = df.iloc[1:, 3]
