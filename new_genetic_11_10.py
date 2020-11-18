# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 15:52:09 2020

@author: StarkPC
"""

### Reference website for genetic algo using roulette wheel selection:
### https://github.com/data-cat-ghub/datacat_blog

import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import mpld3
# from gurobipy import Model, GRB, quicksum
import collections
# import parameters
import time

start_time = time.time()
random.seed(210)

class Simulation:

    def Module_Initialization(self):
        # reading files
        self.df_demand = pd.read_csv("demand.csv")
        self.df_vehicle = pd.read_csv("vehicle_params.csv")
        
        # input vehicle parameters
        self.v_type_dict = {1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:2, 8:2, 9:2, 10:2, 11:3, 12:3, 13:3}
        self.v_cost_dict = {1: 100, 2: 150, 3: 200}
        self.v_max_load_dict = {1: 2500, 2: 3000, 3: 3500} 
        self.v_unit_fuel_cost_dict = {1: 0.2, 2: 0.25, 3: 0.3}
        self.v_c_energy_dict = {1: 4.5, 2: 5, 3: 6}
        self.v_o_energy_dict = {1: 5, 2: 5.5, 3: 6.5}
        self.v_carbon_coeff_dict = {1: 2.75, 2: 2.8, 3: 2.9}
        
        # input demand parameters
        self.d_x_dict = dict(zip(list(self.df_demand["num"]), list(self.df_demand["x_axis_m"])))
        self.d_y_dict = dict(zip(list(self.df_demand["num"]), list(self.df_demand["y_axis_m"])))
        self.d_demand_dict = dict(zip(list(self.df_demand["num"]), list(self.df_demand["demand_per_t"])))
        self.d_e_time_dict = dict(zip(list(self.df_demand["num"]), list(self.df_demand["early_time_hour"])))
        self.d_l_time_dict = dict(zip(list(self.df_demand["num"]), list(self.df_demand["late_time_hour"])))
        self.d_service_time_dict = dict(zip(list(self.df_demand["num"]), list(self.df_demand["service_time_hour"])))
        
        # cumulative travel time to each point
        self.w_time_dict = dict(zip(list(self.df_demand["num"]), [0]*21))

    ###############################################################################
    def Module_Vehicle_Quantity(self):
        list1 = []
        for val in self.vehicle_order_list:
            list1.append(self.v_max_load_dict[self.v_type_dict[val]])
        return list1
    
    ###############################################################################
    def Module_Demand_Quantity(self):
        list1 = []
        for val in self.demand_order_list:
            list1.append(self.d_demand_dict[val])
        return list1
        
    ###############################################################################        
    def Module_Assignment(self, gene_sequence):
        
        self.vehicle_order_list = gene_sequence[1]
        self.vehicle_max_load_list = self.Module_Vehicle_Quantity()
        # eg of vehicle_order_list: [9, 8, 4, 7, 13, 5, 11, 3, 10, 6, 1, 12, 2]
        # eg of vehicle_max_load_list: [3, 3, 3.5, 2.5, 3, 3.5, 2.5, 2.5, 2.5, 3.5, 2.5, 3, 2.5]
    
        self.demand_order_list = gene_sequence[0]
        self.demand_quant_list = self.Module_Demand_Quantity()
            
    ###############################################################################
    def Module_Route_Planner(self):
        
        l1 = []; l2 = []
        l_overall = []; l2_overall = []
        cum_sum = 0; car = 0
        
        vqt_copy = []
        vqt_copy = self.vehicle_max_load_list + [3]*10
        for i in range(len(self.demand_quant_list)):
            dem = self.demand_quant_list[i]
            cum_sum += dem
            if (cum_sum > vqt_copy[car]):
                l_overall.append(l1)
                l2_overall.append(l2)
                l1 = []; l2 = []
                cum_sum = dem
                l1.append(dem)
                l2.append(self.demand_order_list[i])
                car += 1         
            else:
                l1.append(dem)
                l2.append(self.demand_order_list[i])
        l_overall.append(l1)
        l2_overall.append(l2)
        self.route = l2_overall
        
    ###############################################################################
    def Module_Best_Route_Planner(self, gene_sequence):
        l1 = []; l2 = []
        l_overall = []; l2_overall = []
        cum_sum = 0; car = 0
        
        vqt_copy = []
        vqt_copy = self.vehicle_max_load_list + [3]*10
        for i in range(len(self.demand_quant_list)):
            dem = self.demand_quant_list[i]
            cum_sum += dem
            if (cum_sum > vqt_copy[car]):
                l_overall.append(l1)
                l2_overall.append(l2)
                l1 = []; l2 = []
                cum_sum = dem
                l1.append(dem)
                l2.append(self.demand_order_list[i])
                car += 1         
            else:
                l1.append(dem)
                l2.append(self.demand_order_list[i])
        l_overall.append(l1)
        l2_overall.append(l2)
        self.route = l2_overall
        return self.route
        
    ###############################################################################
    def Module_Two_Point_Distance(self, pt1, pt2):
        x1 = self.d_x_dict[pt1]
        x2 = self.d_x_dict[pt2]
        y1 = self.d_y_dict[pt1]
        y2 = self.d_y_dict[pt2]
        
        distance = np.sqrt(((x2-x1)**2) + ((y2-y1)**2))
        return distance
    
    ###############################################################################
    def Module_Two_Point_Time(self, pt1, pt2):
        time = self.Module_Two_Point_Distance(pt1, pt2)/50
        return time
    
    ###############################################################################
    def Module_Distance_Travelled(self):
        route_appended = self.route
        for i in range(len(route_appended)):
            route_appended[i].insert(0,0)
            route_appended[i].append(0)
        
        distance_travelled_list = []    
        for i in range(len(route_appended)):
            dist_travelled = 0
            for j in range(len(route_appended[i])-1):
                dist_travelled += self.Module_Two_Point_Distance(route_appended[i][j], route_appended[i][j+1])
            distance_travelled_list.append(dist_travelled) 
            
        return distance_travelled_list
    
    ###############################################################################
    def Module_Fixed_Cost(self):
        length = len(self.route)
        cost = 0
        if (length <= 13):
            for val in self.vehicle_order_list[0:length]:
                cost += self.v_cost_dict[self.v_type_dict[val]]
        else:
            cost = 1000000
        return cost
    
    ###############################################################################
    def Module_Fuel_Cost(self):
        self.distance_travelled_list = self.Module_Distance_Travelled()
        cost = 0

        if (len(self.route) <= 13):
            for i in range(len(self.route)):
                v_number = self.vehicle_order_list[i]
                v_type = self.v_type_dict[v_number]
                unit_fuel_cons = self.v_unit_fuel_cost_dict[v_type]
                dist_travelled = self.distance_travelled_list[i]
                cost += unit_fuel_cons * dist_travelled
            
            fuel_cost =  6.7 * cost
        else:
            fuel_cost = 1000000
        return fuel_cost
    
    ###############################################################################
    def Module_Carbon_Cost(self):
        #distance_travelled_list = Module_Distance_Travelled(route, vehicle_order_list)
        
        cost = 0
        if (len(self.route) <= 13):
            for i in range(len(self.distance_travelled_list)):
                unit_fuel_cons = self.v_unit_fuel_cost_dict[self.v_type_dict[self.vehicle_order_list[i]]]
                carb_coeff = self.v_carbon_coeff_dict[self.v_type_dict[self.vehicle_order_list[i]]]
                cost += carb_coeff * unit_fuel_cons * self.distance_travelled_list[i]
            
            fuel_cost =  0.5 * cost
        else:
            fuel_cost = 1000000
        return fuel_cost
    ###############################################################################
    def Module_Penalty_Cost(self):
        cum_time = 0
        if (len(self.route) <= 13):
            for pt in self.df_demand["num"]:
                if (pt == 0):
                    continue
                late_time = max((self.w_time_dict[pt] - self.d_l_time_dict[pt]), 0)
                cum_time += late_time
            pen_cost = 2 * cum_time * 60   # hour to minute conversion
        else:
            pen_cost = 1000000
        
        return (pen_cost)
            
    ###############################################################################
    def Module_Damage_Cost(self):
        cum_cost = 0
        if (len(self.route) <= 13):
            for pt in self.df_demand["num"]:
                # skip the inventory location
                if (pt == 0):
                    continue
                cost = (self.d_demand_dict[pt] * (math.exp(1/60 * self.w_time_dict[pt]) - 1))
                cum_cost += cost
            damage_cost = 1 * cum_cost
        else: 
            damage_cost = 1000000
        
        return damage_cost
        
    ###############################################################################
    def Module_Refrigeration_Cost(self):
        self.travel_times_list = [x/50 for x in self.distance_travelled_list]
        
        if (len(self.route) <= 13): 
            cost_travel = 0
            cost_service = 0
            cost_waiting = 0
            # travel time ref costs
            for i in range(len(self.route)):
                travel_time = self.travel_times_list[i]
                v_type = self.v_type_dict[self.vehicle_order_list[i]] 
                v_c_energy = self.v_c_energy_dict[v_type]
                cost_travel += (v_c_energy * travel_time)
            # service time ref costs
            for i in range(len(self.route)):
                route_single = self.route[i]
                v_type = self.v_type_dict[self.vehicle_order_list[i]] 
                v_o_energy = self.v_o_energy_dict[v_type]
                for pt in route_single:
                    service_time = self.d_service_time_dict[pt]
                    cost_service += (v_o_energy * service_time)
            # waiting time ref costs
            for i in range(len(self.route)):
                route_single = self.route[i]
                v_type = self.v_type_dict[self.vehicle_order_list[i]] 
                v_c_energy = self.v_c_energy_dict[v_type]
                
                for pt in route_single:
                    wait_time = max(self.d_e_time_dict[pt] - self.w_time_dict[pt], 0)
                    cost_waiting += v_c_energy * wait_time
                
            ref_cost = 3 * (cost_travel + cost_service + cost_waiting)
        else:
            ref_cost = 1000000
            
        return ref_cost
    
    ###############################################################################
    def Module_Cum_Time_Calculate(self):
        route_appended = self.route
        for i in range(len(route_appended)):
            route_appended[i].insert(0,0)
            
        for i in range(len(route_appended)):
            route_single = route_appended[i]
            
            for i in range(len(route_single)):
                if (route_single[i] == 0):
                    continue
                else:
                    prev_pt = route_single[i-1]
                    current_pt = route_single[i]
                    if (prev_pt == 0):
                        travel_time_recent = self.Module_Two_Point_Time(prev_pt, current_pt)
                        self.w_time_dict[current_pt] = travel_time_recent
                    else:
                        travel_time_recent = self.Module_Two_Point_Time(prev_pt, current_pt)
                        service_time_previous = self.d_service_time_dict[prev_pt]
                        waiting_time_previous = max(self.d_e_time_dict[prev_pt] - self.w_time_dict[prev_pt], 0)
                        self.w_time_dict[current_pt] = self.w_time_dict[prev_pt] + waiting_time_previous + \
                                                        service_time_previous + travel_time_recent
        
    ###############################################################################
    def Module_Calculate_Cost(self, gene_sequence):
        
        self.Module_Initialization()
        # Assign modules
        self.Module_Assignment(gene_sequence)
        self.Module_Route_Planner()
        
        
        if (len(self.route) <= 13):
            self.Module_Cum_Time_Calculate()
            
            fixed_cost  = self.Module_Fixed_Cost()
            fuel_cost   = self.Module_Fuel_Cost()
            carbon_cost = self.Module_Carbon_Cost()
            pen_cost    = self.Module_Penalty_Cost()
            damage_cost = self.Module_Damage_Cost()
            ref_cost    = self.Module_Refrigeration_Cost()
            
            total_cost = fixed_cost + fuel_cost + carbon_cost + pen_cost + damage_cost + ref_cost
            return total_cost,fixed_cost,fuel_cost,carbon_cost,pen_cost,damage_cost,ref_cost
    
        else:
            return 6000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000

###############################################################################
# main part of the code
# ========================================================
# ============ DEFINING MAIN FUNCTIONS (OUTSIDE Sim CLASS) ============
# ========================================================
###############################################################################
# initial decision variables
sim = Simulation()
decision_variables = [list(range(1,21)), list(range(1,14))]
n_children = 40
colors_dict = {1:'b', 2:'g', 3:'r'} 
shapes_list = ["o", "v", "s", "^", "<", "P"]
def initialize():
    pop_bag = [[[], []], [[], []], [[], []], [[], []], [[], []], 
               [[], []], [[], []], [[], []], [[], []], [[], []],
               [[], []], [[], []], [[], []], [[], []], [[], []], 
               [[], []], [[], []], [[], []], [[], []], [[], []],
               [[], []], [[], []], [[], []], [[], []], [[], []], 
               [[], []], [[], []], [[], []], [[], []], [[], []],
               [[], []], [[], []], [[], []], [[], []], [[], []], 
               [[], []], [[], []], [[], []], [[], []], [[], []]]
    # demand
    for i in range(n_children):
        for j in range(len(decision_variables)):
            rnd_sol = decision_variables[j].copy()
            random.shuffle(rnd_sol)
            pop_bag[i][j] = (rnd_sol)
    return np.array(pop_bag)

###############################################################################
def fitness_function(solution):
    all_cost = sim.Module_Calculate_Cost(solution)
    total_cost = all_cost[0]
    return total_cost

###############################################################################
def eval_fit_population(pop_bag):
    result = {}
    fit_vals_lst = []
    solutions = []
    for solution in pop_bag:
        fit_vals_lst.append(fitness_function(solution))
        solutions.append(solution)
    result["fit_vals"] = fit_vals_lst
    min_wgh = [np.max(list(result["fit_vals"]))-i for i in list(result["fit_vals"])]
    result["fit_wgh"]  = [val/sum(min_wgh) for val in min_wgh]
    result["solution"] = np.array(solutions)
    return result

###############################################################################
def pickOne(pop_bag):
    fit_bag_evals = eval_fit_population(pop_bag)
    a=True
    while a:
        rnIndex = random.randint(0, len(pop_bag)-1)    # length of pop_bag = 10
        rnPick  = fit_bag_evals["fit_wgh"][rnIndex]
        r = random.random()   # between 0 and 1
        if  r <= rnPick:
            pickedSol = fit_bag_evals["solution"][rnIndex]
            a = False
    return pickedSol

###############################################################################
def crossover(solA, solB):
    
    d_solA = solA[0]
    d_solB = solB[0]
    v_solA = solA[1]
    v_solB = solA[1]
    #######################
    n1 = len(d_solA)
    d_child = [np.nan for i in range(n1)]
    d_num_els = np.ceil(n1*(random.randint(10,90)/100))
    d_str_pnt = random.randint(0, n1-2)
    d_end_pnt = n1 if int(d_str_pnt+d_num_els) > n1 else int(d_str_pnt+d_num_els)
    d_blockA = list(d_solA[d_str_pnt:d_end_pnt])
    d_child[d_str_pnt:d_end_pnt] = d_blockA
    for i in range(n1):
        if list(d_blockA).count(d_solB[i]) == 0:
            for j in range(n1):
                if np.isnan(d_child[j]):
                    d_child[j] = d_solB[i]
                    break
    #######################
    n2 = len(v_solA)
    v_child = [np.nan for i in range(n2)]
    v_num_els = np.ceil(n2*(random.randint(10,90)/100))
    v_str_pnt = random.randint(0, n2-2)
    v_end_pnt = n2 if int(v_str_pnt+v_num_els) > n2 else int(v_str_pnt+v_num_els)
    v_blockA = list(v_solA[v_str_pnt:v_end_pnt])
    v_child[v_str_pnt:v_end_pnt] = v_blockA
    for i in range(n2):
        if list(v_blockA).count(v_solB[i]) == 0:
            for j in range(n2):
                if np.isnan(v_child[j]):
                    v_child[j] = v_solB[i]
                    break    
    child = [d_child, v_child]
    return child

###############################################################################
def mutation(sol):
    n1 = len(sol[0])
    d_pos_1 = random.randint(0,n1-1)
    d_pos_2 = random.randint(0,n1-1)
    d_result = swap(sol[0], d_pos_1, d_pos_2)
    #######################
    n2 = len(sol[1])
    v_pos_1 = random.randint(0,n2-1)
    v_pos_2 = random.randint(0,n2-1)
    v_result = swap(sol[1], v_pos_1, v_pos_2)
    result = [d_result, v_result]
    return result

###############################################################################
def swap(sol, posA, posB):
    d_result = sol.copy()
    elA = sol[posA]
    elB = sol[posB]
    d_result[posA] = elB
    d_result[posB] = elA
    return d_result

###############################################################################
def Closed_Route(best_route_1):
    list_2 = []
    for i in range(len(best_route_1)):
        single_route = best_route_1[i]
        single_route.insert(0,0)
        single_route.append(0)
        list_2.append(single_route)

    print (list_2)
    return (list_2)

###############################################################################                  
def Plotter(best_route_1):
    k = 0
    for i in range(len(best_route_1)):
        single_route = best_route_1[i]
        for j in range(len(single_route)-1):
            cur_pt = single_route[j]
            nxt_pt = single_route[j+1]
            point1 = [sim.d_x_dict[cur_pt], sim.d_y_dict[cur_pt]]
            point2 = [sim.d_x_dict[nxt_pt], sim.d_y_dict[nxt_pt]]
            x_values = [point1[0], point2[0]]
            y_values = [point1[1], point2[1]]
            
            v_number = sim.vehicle_order_list[i]
            v_type = sim.v_type_dict[v_number]
            color_plot = colors_dict[v_type]
            plt.plot(x_values, y_values, color=color_plot, linestyle='dashed', marker=shapes_list[k%6])
            
        k+=1
   
###############################################################################
# ========================================================
# ============ START THE EVOLUTIONARY PROCESS ============
# ========================================================

# Create the initial population bag
pop_bag  = initialize()
overall_best_fit = []
# Iterate over all generations
for g in range(10):
    print ("generation = ", g)
    # Calculate the fitness of elements in population bag
    pop_bag_fit = eval_fit_population(pop_bag)
    
    # Best individual in the current population bag
    best_fit = np.min(pop_bag_fit["fit_vals"])
    best_fit_index = pop_bag_fit["fit_vals"].index(best_fit)
    best_solution  = pop_bag_fit["solution"][best_fit_index]
    
    # Check if we have a new best
    if g == 0:
        best_fit_global      = best_fit
        best_solution_global = best_solution
    else:
        if best_fit <= best_fit_global:
            best_fit_global      = best_fit
            best_solution_global = best_solution
    overall_best_fit.append(best_fit_global)

    print ("best cost = ", best_fit_global)
    # Create the new population bag
    new_pop_bag = []

    for i in range(n_children):
        # Pick 2 parents from the bag
        pA = pickOne(pop_bag)
        pB = pickOne(pop_bag) 
        new_element = pA
        
        # Crossover the parents
        if random.random() <= 0.87:
            new_element = crossover(pA, pB)
        # Mutate the child
        if random.random() <= 0.7:
            new_element = mutation(new_element) 
        # Append the child to the bag
        new_pop_bag.append(new_element)

    # Set the new bag as the population bag
    pop_bag = np.array(new_pop_bag)

# Best fitness and solution
print(f"Best Fitness: {best_fit_global}")
print(f"Best Solution: {best_solution_global}")

# Post processing
sim.Module_Assignment(best_solution_global)
best_route = sim.Module_Best_Route_Planner(best_solution_global)
print ("best route = ", best_route)
best_route_1 = Closed_Route(best_route)
print("best route = ", best_route_1)
Plotter(best_route_1)

print("--- %s seconds ---" % (time.time() - start_time))
