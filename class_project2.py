import os
import numpy as np
import re
import scipy 


def split_equation(equation):
  # Regular expression to split and avoid numerics as unique elements
  reactants = re.findall(r"([A-Za-z]+)", equation.split("->")[0] 
                           if "->" in equation else equation.split("=")[0])
  products = re.findall(r"([A-Za-z]+)", equation.split("->")[1] 
                          if "->" in equation else equation.split("=")[1])
  return reactants , products
    
def y_split_equation(equation):
  reactants = (equation.split("->")[0] 
               if "->" in equation else equation.split("=")[0]).split("+")
  products = (equation.split("->")[1] 
              if "->" in equation else equation.split("=")[1]).split("+") 
  return reactants , products

def Z_split_equation(equation):
  reactants = re.findall(r"[HO]", equation.split("->")[0] 
                           if "->" in equation else equation.split("=")[0])
  products = re.findall(r"[HO]", equation.split("->")[1] 
                          if "->" in equation else equation.split("=")[1])
  return reactants , products

eq = 'H2 + O2 = CO2 + H2O'
equation = 'H2 + O2 = H2O'

#print(split_equation(eq))
reacts, prods = Z_split_equation(eq)
reacts = [x.strip() for x in reacts]
prods = [x.strip() for x in prods]
u_reacts = list(set("".join( reacts + prods)))
#print(u_reacts)
#list_of_strings = ['a','a','c','d','e','c','f']
#uniques = set(list_of_strings)
#print(uniques)

#coefficient matrix
matrix_array = np.zeros((len(u_reacts),len(reacts + prods))) #dimension of the matrix (())
#populating the matrix with the respective items
for i,u in enumerate(u_reacts): #lists the elements and their indexes
    #print(i,u)
    for j,cpd in enumerate(reacts + prods):
        #print(j,cpd)
        num_of_elems = cpd.count(u) #count the number of times u appears in k
        #print(num_of_elems)
        matrix_array[i, j] = num_of_elems #populate the matrix with the indexes of the u & cpds 
print(matrix_array) #rows = u columns = cpd if element in the row appears in the column or not
      
#soln vector
h_vector = np.zeros(len(reacts) + len(prods))
h_vector[len(reacts):] = 1 #sets all elements of prods to 1 ie there's prods

#solving system of linear eqs
x = np.linalg.lstsq(matrix_array.T, h_vector, rcond=None)[0]
print(x)


# Print the balanced equation
balanced_equation = ''
for j, u in enumerate(u_reacts):
  coeff = x[j]
  if coeff == 1:
    balanced_equation += u
  else:
    balanced_equation += '{:.0f}{}'.format(coeff, u)
  if j < len(u_reacts) - 1:
    balanced_equation += ' + '
print(balanced_equation)



