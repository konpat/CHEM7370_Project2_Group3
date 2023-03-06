import os
import sys
import numpy as np
from scipy import linalg
from fractions import Fraction

######## This script is used to balance molecular or ionic chemical reactions. 
######## At this point, it is not too effective for redox reactions that contain e- in the equation. 

######## execute the script with balance('reactions.txt')

############################################


### count the number of atoms in a molecule. Taken from the course material CHEM7370 Class 7 Class7solved.ipynb
def count_atoms(molecule):
    atoms = {}
    current_atom = ''
    current_number = ''
    for x in molecule:
        if x.isupper() or x == '+' or x == '-':
# x is a capital letter - start of a new atom
            if current_atom != '':
                if current_number == '': 
                    if current_atom in atoms.keys():
                        atoms[current_atom] += 1
                    else:
                        atoms[current_atom] = 1
                else:
                    if current_atom in atoms.keys():
                        atoms[current_atom] += int(current_number)
                    else:
                        atoms[current_atom] = int(current_number)
                current_atom = ''
                current_number = ''
            current_atom += x
        elif x.islower():
# x is a lowercase letter - second part of an atom
            current_atom += x
        else:
# x is a number - part of the coefficient
            current_number += x
        
    #print(current_atom,current_number)
    if current_atom != '':
        if current_number == '':
            if current_atom in atoms.keys():
                atoms[current_atom] += 1
            else:
                atoms[current_atom] = 1
        else:
            if current_atom in atoms.keys():
                atoms[current_atom] += int(current_number)
            else:
                atoms[current_atom] = int(current_number)
        
    return atoms


### construct the matrix ingredients of the reaction system
def arr(lib, libname):
    
    arr_all = []

    for m in range(len(lib)):
#        print(libr[m])
        arr = []
        uplibname = libname.copy()
             
        for k, v in lib[m].items():
            if k == '+':
                uplibname['e'] -= v
                
            elif k == '-':
                uplibname['e'] += v
            
            else: 
                uplibname[k] += v
                
        arr = list(uplibname.values())
        arr_all.append(arr)
        
#        print(arr_all)       
#        print(uplibname)

    return arr_all


### convert the raw coefficients into the fraction format and find the smallest denominator
def fraction(array):
    
    uparray = []
    
    fac = []
    
    for i in range(len(array)):
        t = ()
        f = str(Fraction(array[i]).limit_denominator(1000))
        if '/' not in f:
            t = (float(f), 0.0)
        else:         
            t = (float(f.split('/')[0]),float(f.split('/')[1]))
        
        uparray.append(t)
        
    for i in range(len(uparray)):        
        a = uparray[i][1]
        
        if a > 0.0:
            fac.append(a)
    
    if fac != []:
        minfac = min(fac)
    else:
        minfac = 1.0
            
    return uparray, minfac


### calculate the final coefficients for the balanced reaction
def coeff_arr(array):
    
    coe = []
    ficoe = []
    
    coe = np.insert(array, 0, [1])
    
    arr_frac, factor = fraction(coe)
    #print(arr_frac, factor)
    
    arr_frac_1 = factor * coe   
    #print(arr_frac_1)
    
    arr_frac_fi, factor_1 = fraction(arr_frac_1)  
    #print(arr_frac_fi, factor_1)
      
    for i in range(len(arr_frac_fi)):
        a = int(abs(arr_frac_fi[i][0]))
        ficoe.append(a)
            
    return ficoe


### calculate the remain charge of the reactants and products
def remain_charge(plus, minus):

    re = ()
    
    if plus > minus:
        re = ('+', plus - minus)
    elif minus > plus:
        re = ('-', minus - plus)
    else:
        re = (0.0, 0.0)
        
    return re


### check the status of the reaction what is balanced so far
def check(coeff, libr, libp, libname):
    
    t = ''
    
    coeff_reacts = coeff[:len(libr)]
    coeff_prods = coeff[len(libr):]
    
    charge = {'+': 0, '-': 0}
    
    if 'e' in libname.keys():
        libname.pop('e')    
   # print(libname)
       
    libname_r = libname.copy()
    libname_r.update(charge)
    
    libname_p = libname_r.copy()
    
  #  print(libname_r)
  #  print(libname_p)
    
    for i in range(len(libr)):
        for m, n in libr[i].items():             
            libname_r[m] += n * coeff_reacts[i]
    
    for i in range(len(libp)):
        for m, n in libp[i].items():             
            libname_p[m] += n * coeff_prods[i]
                
   # print(libname_r)
   # print(libname_p)
    
    residue_charge_r = remain_charge(libname_r['+'], libname_r['-'])    
    residue_charge_p = remain_charge(libname_p['+'], libname_p['-'])
    
   # print(residue_charge_r)
   # print(residue_charge_p)
    
    libname_r.pop('+')
    libname_r.pop('-')
    
    libname_p.pop('+')
    libname_p.pop('-')
    
    for k, v in libname_r.items():
        for i, j in libname_p.items():
            if k == i and v == j:
                if residue_charge_r == residue_charge_p:
                    t = 'molecule and charges are balanced!'
                else:
                    t = 'charges are not balanced!'
            else:
                
                if residue_charge_r == residue_charge_p:
                    t = 'molecule are not balanced!'
                else:
                    t = 'molecule and charges are not balanced!'
    
    return t


### prepare the answer for the balanced reaction
def final(coeff, reacts, prods):
    
    r = ''    
    rp = ''    
    rpf = ''
    
    coeff_reacts = coeff[:len(reacts)]
    coeff_prods = coeff[len(reacts):]
    
    for i in range(len(reacts)):
                
        if coeff_reacts[i] != 1:
            r += str(coeff_reacts[i]) + reacts[i] + ' ' + '+' + ' '
        else:
            r += reacts[i] + ' ' + '+' + ' '
    
    rp = r[::-1].replace('+','>-',1)[::-1]

    for i in range(len(prods)):
        
        if coeff_prods[i] != 1:
            rp += str(coeff_prods[i]) + prods[i] + ' ' + '+' + ' '
        else:
            rp += prods[i] + ' ' + '+' + ' '
            
    rpf = rp[::-1].replace('+','',1)[::-1]
    
    return rpf
            

    
def balance(myfile):
        
    with open(myfile) as f:
        lines = f.readlines()
        
        for m in range(len(lines)):
            
            reaction = ''
    
            libr = {}
            libp = {}
            
            libname = {}

            reaction = str(lines[m])            
        #    print(reaction)
            
            react_all = reaction.split(" -> ")[0]
            prod_all = reaction.split(" -> ")[1]
            reacts = react_all.split(" + ")                      ## reactant molecules
            prods = prod_all.split(" + ")
                
            prods = [prods[i].replace('\n', '') for i in range(len(prods))]            ## product molecules

        #    print(reacts)
        #    print(prods)

            for i in range(len(reacts)):
                libr[i] = count_atoms(reacts[i])                 ## library for each molecule 
                
                for j in libr[i].keys():                       
                    if j not in libname.keys():
                        if j != '+' and j != '-':
                            libname[j] = 0 
            
            for i in range(len(libr.keys())):
                for k, v in libr[i].items():
                    if k == '+' or k == '-':
                        libname['e'] = 0                         
            
       #     print(libr)
       #     print(libname)
            
            for i in range(len(prods)):
                libp[i] = count_atoms(prods[i])
       #         print(libp[i])
                
            for i in range(len(libp.keys())):
                for k, v in libp[i].items():
                    if k == '+' or k == '-':
                        libname['e'] = 0                        ## collect unique atoms and charge present in the reaction
            
       #     print(libp)
       #     print(libname)
            
            arrr = arr(libr, libname)            
            arrp = arr(libp, libname)                        
       #     print(arrr)
       #     print(arrp)
            
            arr_sum = arrr + arrp                                          
       #     print(np.array(arr_sum))                    
       #     print(np.array(arr_sum).T)
            
            mat = np.array(arr_sum[1:]).T                      ## matrix of the coefficients system
       #     print(mat)
            
            vec = np.array(arr_sum[0])
       #     print(vec)
            
            ans = np.linalg.lstsq(mat, vec, rcond=None)[0]     ## coefficients of the reaction 
       #     print(ans)
            
            ficoe = coeff_arr(ans)                             ## modified coefficients 
       #     print(ficoe)

            status = check(ficoe, libr, libp, libname)         ## check the status of the balanced reaction
       #     print(status)
            
            if status == 'molecule and charges are balanced!':
                balanced_reaction = final(ficoe, reacts, prods)          ## the balanced reaction
                
            else:               
                balanced_reaction = status
            
            r = open('balanced.txt', 'a+')
            
            r.write("%s \n" % balanced_reaction)
            
            r.close()
            
    return balanced_reaction
           



