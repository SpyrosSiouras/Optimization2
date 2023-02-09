import numpy as np
import pandas as pd
from math import log10
import csv



def rosenbrock(x):
    global function_calls
    function_calls +=1
    return (1-x[0])**2 + 100*(x[1]-x[0]**2)**2

def predator_prey(x):
    global function_calls
    function_calls +=1
    p = data[0].to_list()
    q = data[1].to_list()

    p_pred = [1]
    q_pred = [1]

    error_p = [0]
    error_q = [0]
    r = x[0]
    K = x[1]
    s = x[2]
    u = x[3]
    v = x[4]
    max_error = 0
    currentP = 1
    currentQ = 1
    for i in range(1, 101):
        deltaP = (r * (1 - currentP /  K) -  s*currentQ) * currentP
        deltaQ = (-u + v*currentP) * currentQ
        # print(f"{(1 - currentP /  K)=}")
        # print(f"{deltaQ=}")
        currentP = currentP + deltaP
        currentQ = currentQ + deltaQ

        p_pred.append(currentP)
        q_pred.append(currentQ)

        errorP = abs(p[i] - currentP) / p[i]
        errorQ = abs(q[i] - currentQ) / q[i]

        error_p.append(errorP)
        error_q.append(errorQ)

        if errorP > max_error:

            max_error = errorP

        if errorQ > max_error:

            max_error = errorQ

    return max_error

def nelder_mead(func, x0, bounds,alpha=1, gamma=2, rho=-0.5, sigma=0.5, epsilon=1e-6, max_iter=4000):
    n = len(x0)
    simplex = np.zeros((n+1, n))
    simplex[0] = x0
    for i in range(1, n+1):
        simplex[i] = x0 + alpha * np.eye(n)[i-1]
    
    for i in range(1,max_iter):
        print(f"Iteration: {i}, best so far = {simplex[0]}, functioncalls: {function_calls}")
        simplex = np.array(sorted(simplex, key=lambda x: func(x)))
        x_bar = np.mean(simplex[:-1], axis=0)

        x_r = x_bar + gamma * (x_bar - simplex[-1])

        f_r = func(x_r)


        if f_r < func(simplex[0]):
            x_e = x_bar + rho * (x_r - x_bar)
            f_e = func(x_e)
            if f_e < f_r:
                simplex[-1] = x_e
            else:
                simplex[-1] = x_r
        elif func(simplex[-2]) < f_r < func(simplex[-1]):
            x_c = x_bar + rho * (x_r - x_bar)
            f_c = func(x_c)
            if f_c < func(simplex[-1]):
                simplex[-1] = x_c
            else:
                simplex[-1] = x_bar + sigma * (simplex[-1] - x_bar)
                simplex[-2] = x_bar + sigma * (simplex[-2] - x_bar)
        else:
            x_c = x_bar - rho * (x_bar - simplex[-1])
            f_c = func(x_c)
            if f_c < func(simplex[-1]):
                simplex[-1] = x_c
            else:
                simplex[-1] = x_bar + sigma * (simplex[-1] - x_bar)
                simplex[-2] = x_bar + sigma * (simplex[-2] - x_bar)
        

        if func(simplex[0]) < epsilon or function_calls > 1_000_000:
            print(f"Reached wanted accuracy!")
            print()
            break
        
        if i % 500 == 0:
            
            coef = (log10(i+2))
            for j in range(3): 
                for k in range(len(x0) - 1):

                    simplex[j, k] += np.random.uniform(low= (-1  / coef) , high= (1  / coef) )
                
        
        if i > 3000: 
           
            x0 = np.random.uniform(low=0.1, high=3,size= len(x0))
            return (nelder_mead(func, x0, bounds, 0.035, max_iter = 3100))

        for j in range(3):
            for k in range(len(x0)):
                simplex[j,k] =  max(bounds[k][0], min(simplex[j, k], bounds[k][1]))
    
    return simplex[0]


def run(func, x_size, accuracy, max_iterations):
    global function_calls
    x0 = np.random.uniform(low=0.1, high=3,size= x_size)
    print(x0)
    # x0 = np.array([1.01882, 1.0166, 0.513933, 0.684169, 1.55421])
    bounds = [(1.18, 1.37),  (0.93, 1.08), (0.42, 0.60), (0.6, 0.8),(1.52, 1.78)]
    x_min = nelder_mead(func, x0, bounds,epsilon=accuracy, max_iter= max_iterations)

    print(f"{x_min=}")
    print(f"{func(x_min)}")
    print(f"{function_calls=}")

    return [func(x_min), x_min[0], x_min[1],x_min[2],x_min[3],x_min[4], function_calls]


if __name__ == "__main__":
    global function_calls
    
    function_calls = 0
    data = pd.read_csv("data.csv", header=None)
    # run(rosenbrock, 2, 1e-6, 100)


    for _ in range(10):

        results = run(predator_prey, x_size = 5, accuracy = 0.035, max_iterations = 6000)


        with open("nelderMead.csv", 'a', encoding='UTF8', newline='') as f:

            writer = csv.writer(f)

            writer.writerow(results)

        function_calls = 0