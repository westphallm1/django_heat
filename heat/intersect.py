import solve_eq
import numpy as np

#constants
R = 8.314
Cv = 1.5 * R
Cp = 2.5 * R
gamma = Cp/Cv

def A_str(Vo,Po):
    """
    Returns pressure as a function of volume for an adiabat at the given P & V
    P = Po * Vo**gamma / V**gamma
    """
    return "%f / x**%f"%(Po*Vo**gamma,gamma)
    
def B_str(Vo,Po):
    """
    Returns pressure as a function of volume for an isobar at the given P & V
    P = Po
    """
    return "%s"%Po
    
#the volume and pressure at an isochor are known, handled in find_intersect
    
def T_str(Vo,Po):
    """
    Returns pressure as a function of volume for an isobar at the given P & V
    P = Po*Vo/V
    """
    return "%f / x"%(Po*Vo)
    
def eqn_str(Vo,Po,eqn_type):
    """calls the appropriate X_str function"""
    if eqn_type == "A":
        return A_str(Vo,Po)
    if eqn_type == "B":
        return B_str(Vo,Po)
    if eqn_type == "C":
        return None
    if eqn_type == "T":
        return T_str(Vo,Po)
        
        
def find_intersect(V0,P0,Vf,Pf,type0,typef):
    """Find the volume at which pressure the thermodynamic processes 
    branching from (Vo, Po) and (Vf, Pf) intersect
    V and P are the pressure and volume at each point
    type0 and typef are the process occuring:
      "A": adiabatic
      "B": isobaric
      "C": isochoric
      "T": isothermal 
    returns the volume and pressure of intersect
    """
    
    #get the equation strings for each process
    eqn0 = eqn_str(V0,P0,type0)
    eqnf = eqn_str(Vf,Pf,typef)
    
    #Find the volume and pressure at the intersect of the curves
    Vm = None
    Pm = None
    
    #if an isochoric process is occuring, the volume of intersection
    #is already known
    if type0 == "C":
        Vm = V0
        Pm = solve_eq.valfstr(eqnf,V0)
        #check that Pm lies with P0 and Pf
        if Pm < P0 or Pm > Pf:
            Pm = Vm = None
    elif typef == "C":
        Vm = Vf
        Pm = solve_eq.valfstr(eqn0,Vf)
        #check that Pm lies with P0 and Pf
        if Pm < P0 or Pm > Pf:
            Pm = Vm = None
    else:
        #solve the equation
        eqn = ' - '.join((eqnf,eqn0))
        Vm = solve_eq.find_0(eqn,Vf,V0)
        if Vm is None: 
            Pm = None
        else: 
            Pm = solve_eq.valfstr(eqn0,Vm)
    
    return Vm, Pm
