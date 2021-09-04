import pandas as pd
import numpy as np
data_U = pd.read_csv("ieee39节点电压及功率.csv", encoding='utf-8',index_col=0)

U = np.array(data_U["V"])
S = np.array(data_U["S_node"])
#
# 读取导纳矩阵
Y = pd.read_csv("Y_None_None.csv",index_col=0).to_numpy()

# 把读取的字符串转化成复数
for i,line in enumerate(Y) :
    for j,row in enumerate(line):
        Y[i][j]=complex(Y[i][j])

for i,line in enumerate(U) :
    U[i]=complex(line)
# 读取节点数据
bus_data_list =pd.read_csv('bus_data.csv',index_col=0)
bus_data_list.set_index(['Bus number'],inplace=True)
node_list = bus_data_list.index.values
# print(node_list)
# 各节点自导纳

gen_x_dict = {30:complex(0,0.031),31:complex(0,0.0697),32:complex(0,0.0531),33:complex(0,0.0436),34:complex(0,0.132),35:complex(0,0.05),36:complex(0,0.049),37:complex(0,0.057),38:complex(0,0.057),39:complex(0,0.006)}
# gen_x_loc_list = list(gen_x_dict.keys())
gen_num = len(gen_x_dict)
bus_num = len(node_list)
for node_num in node_list:
    bus_data = bus_data_list.loc[node_num]
    # 数组的实际索引是节点数-1
    node_num = node_num-1
    # 修改节点导纳矩阵
    # 修改负荷影响的部分

    # 获取该母线数据
    loadP =  float(bus_data['Load MW'])/100
    loadQ = float(bus_data['Load MVAR'])/100
    loadS = complex(loadP,loadQ)
    my_U = U[node_num]
    # 由于负荷接地，修改自导纳即可
    if loadP != 0 or loadQ != 0:
        Z_ld = loadS * (abs(my_U) ** 2) / (abs(loadS) ** 2)
        Y[node_num][node_num] +=  1/Z_ld
    # Y[node]
    # 修改由电机影响的部分
    if node_num+1 in gen_x_dict:
        xt = gen_x_dict[node_num+1]
        Y[node_num][node_num] += 1 / xt
# # 求节点阻抗矩阵
Y = Y.astype(np.complex)
Y = np.around(Y,decimals=2)
pd.DataFrame(Y).to_csv('fault_y1.csv')
# short_data = pd.DataFrame({'I':I,'S':S,'Sn':Sn})
I = np.zeros([len(node_list)],dtype=complex)
S = np.zeros([len(node_list)],dtype=complex)
Sn = np.zeros([len(node_list)],dtype=complex)
Z = np.linalg.inv(Y)
for node_num in node_list:
    # 数组的实际索引是节点数-1
    node_num = node_num - 1
     # 计算短路电流有名值
    I[node_num] =  U[node_num] /  Z[node_num][node_num]
    S[node_num] = abs(I[node_num])
    Sn[node_num] =  S[node_num] * 100
short_data = pd.DataFrame({'I':abs(I),'S':S,'Sn':Sn})
short_data.to_csv('short_data.csv')
#
# print('计算结束')
# print('abs( I)=')
# print(abs( I))
# print('Sn=')
# print( Sn)
#
# if __name__ == '__main__':
#     my_flaut = Flaut()
