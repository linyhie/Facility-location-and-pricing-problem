import numpy as np

def read_the_data(file_name):
    file1 = open(file_name, "r")          
    Lines = file1.readlines()
    Length = len(Lines)    
    for n in range(Length):
        if n == 0:
            J, I, r = int(((Lines[n].strip()).split())[0]), int(((Lines[n].strip()).split())[1]), int(((Lines[n].strip()).split())[2]) 
            c = np.zeros((I,J))
            b = np.zeros(I)
        elif n < Length - 1:
            for i in range(I):
                c[i,n-1] = int(((Lines[n].strip()).split())[i])           
        else:
            for i in range(I):
                b[i] = int(((Lines[n].strip()).split())[i])
    return(I,J,r,c,b)
####################################################

file_name = "FLPr_100_40_01.txt"   ### J = {40,100}, inst = {01,02,...,10}

I,J,r,c,b = read_the_data(file_name)
####################################################

K = 10              # number of pricing levels

timelimit = 7200

price_bound = 60
### pricing level
p = np.zeros((J,K))
for j in range(J):
    p[j] = (np.arange(0,K)+1)*price_bound/K
    
theta = np.zeros((I,J,K))
for i in range(I):
    for j in range(J):
        theta[i,j] = c[i,j] + p[j]

pi = np.zeros((I,J,K))
for i in range(I):
    pi[i] = b[i] - theta[i]
####################################################
### delta_{ijk} means whether facility j at k level is dominated by outside option, 
### delta_{ijk} = 1 if it is dominated, i.e., i.e., pi[i,j,k] = b[i] - theta[i,j,k] < 0 
delta = pi < 0  ### y_{ijk} can be 1 only if pi_{ijk} > 0

### starting tracking time
import time as TIME
start_time = TIME.time()
####################################################
import gurobipy as grb
model = grb.Model()
model.setParam('OutputFlag', 1)
model.setParam('TimeLimit',timelimit)
model.setParam('DisplayInterval',50)

x = model.addVars(J,K,vtype=grb.GRB.BINARY)
y = model.addVars(I,J,K) 
w = model.addVars(J,vtype=grb.GRB.BINARY)

model.setObjective(grb.quicksum(p[j,k]*y[i,j,k] for i in range(I) for j in range(J) for k in range(K)), 
                                grb.GRB.MAXIMIZE)

model.addConstr(grb.quicksum(w[j] for j in range(J)) == r)

for j in range(J):model.addConstr(grb.quicksum(x[j,k] for k in range(K)) == w[j])

for i in range(I):
    model.addConstr(grb.quicksum(y[i,j,k] for j in range(J) for k in range(K)) <= 1)    
    for j in range(J):
        for k in range(K):
            model.addConstr(y[i,j,k] <= x[j,k])
            model.addConstr(y[i,j,k] <= 1 - delta[i,j,k])  ### Preprossing: y_{ijk} can be 1 only if delta[i,j,k] = 0      

def solve_follower_problem(set_open):
     ### solve for y        
     y_last = np.zeros((I,J,K)) ### best y found by solving the follower problem
     number_facilities = len(set_open)         
     for i in range(I):
         utility_of_open = []
         for o in range(number_facilities):
             utility_of_open.append(pi[i][set_open[o]])
         index = np.argmax(utility_of_open)
         if utility_of_open[index] >= 0:          ### select some facility
             y_last[i][set_open[index]] = 1
     return(y_last)
        
        
def lazy_cut(model, where):    
    if where == grb.GRB.Callback.MIPSOL:
         print("--sep--")
         x_vals = model.cbGetSolution(model._x) 
         y_vals = model.cbGetSolution(model._y)                  
         set_open = [(j,k) for j in range(J) for k in range(K) if x_vals[j,k] > 0.5]
         x_sol = np.zeros((J,K))    ### x solution from the tree      
         y_sol = np.zeros((I,J,K))  ### y solution from the tree
         for j in range(J):
             for k in range(K):
                 x_sol[j,k] = x_vals[j,k]
                 for i in range(I):
                     y_sol[i,j,k] = y_vals[i,j,k]
         ### solve for y        
         y_last = solve_follower_problem(set_open)
         ## add cuts                            
         piy = pi*y_last         
         for i in range(I): 
             if np.sum(pi[i]*y_sol[i]) < np.sum(piy[i]*x_sol): ### check cuts violation
                 model.cbLazy(grb.quicksum(pi[i,j,k]*y[i,j,k] for j in range(J) for k in range(K))
                           >= grb.quicksum(piy[i,j,k]*x[j,k] for j in range(J) for k in range(K)))
                 model._number_of_cuts +=1

model._x = x
model._y = y
model._number_of_cuts = 0
model.Params.lazyConstraints = 1
model.optimize(lazy_cut) 
totla_time = TIME.time() - start_time 
### end tracking time
open_facility = [j for j in range(J) if w[j].x > 0.5]
w_sol = [round(w[j].x) for j in range(J)]
print("********************************************")    
print("Total CPU", round(totla_time,1))    
print("Solver CPU", round(model.Runtime,1))
print("rgap", round(model.MIPgap*100,2))
print("Objevtive", round(model.ObjVal,1))
print("open facilities", open_facility)
x_sol = np.zeros((J,K))
for j in range(J):
    for k in range(K):
        x_sol[j,k] = round(x[j,k].x)
price = np.sum(x_sol*p,axis=1)
print("with price", price[np.nonzero(price)])
print("average price", round(np.average(price[np.nonzero(price)]),1))
print("# cuts",  model._number_of_cuts) 
print("Number of branch and cut nodes is ", round(model.NodeCount))   