from copy import deepcopy
import numpy as np

import json

# a = np.array([[[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]], [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]] ]).T
# b = np.array([2, 3])
a = -3/4 * np.pi
b = np.sin(a)
c = a*b
# print(b)
print(b)
# a = np.array([[[2, 2, 3, 4, 5], [2, 2, 3, 4, 5]], [[2, 2, 3, 4, 5], [2, 2, 3, 4, 5]]])

# print(np.size(a, 1))
 
# for i in range(1, 5):
#     print(i)


# dict1 = {'a':1, 'b': 2}
# dict2 = {'q':3, 'p':4}
# dict_list = []
# dict_list.append(dict1)

# dict_list.append(dict2)



# with open('demo3.json', mode='w') as f:
#     temp = json.dumps(dict_list)
    # json.dump(dict_list, f)
# aaa = []
# llll = []
# arr = np.zeros(1)
# inhaleMass = {}
# for i in range(3):
#     inhaleMass[i] = i
#     if i == 2:
#         inhaleMass[i] = 6
#     arr[0] = i
#     llll.append(deepcopy(arr))
#     aaa.append(deepcopy(inhaleMass))

# print(llll)


# bound = [[123, 324], [5413, 4324]]
 

# net_X = bound[1][0] - bound[0][0]
# net_Y = bound[1][1] - bound[0][1]

# Y = net_X * (1+ 0.1)
# X = net_Y * (1+ 0.1)
# gridSize = 218

# x, y = np.meshgrid(range(round(bound[0][0] - 0.05*net_X), round(bound[1][0] + 0.05*net_X), gridSize), range(round(bound[0][1] - 0.05*net_Y), round(bound[1][1] + 0.05*net_Y), gridSize))
     
# print(np.size(x, 0))
# print(int (X / gridSize) + 1)

# print(np.size(y, 1))
# print(int (Y / gridSize) + 1)


# #print(y)
# print((Y / gridSize))