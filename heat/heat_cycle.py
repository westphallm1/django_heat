import numpy as np
import intersect
import solve_eq
import matplotlib
matplotlib.use("Agg")


class HeatCycle():
    """Stores the endpoints of the processes in a heat engine, as well as the
    types of process, and methods for generating a "valid" cycle
    In all methods in this class:
      "A" refers to an adiabatic process
      "B" refers to an isobaric process
      "C" refers to an isochoric process
      "T" refers to an isothermal process
    """
    def __init__(self,V0,P0,Vf,Pf,nodes = 0, random_default = True):
        """(V0,P0) and (Vf,Pf) are the bottom left and upper right corners of
        the heat cycle, respectively, and nodes is the number of additional
        process endpoints between the two extremes.
        """
        self.V0 = V0
        self.P0 = P0
        self.Vf = Vf
        self.Pf = Pf
        self.nodes = nodes

        #Initialize arrays to store the name of each process in the cycle,
        #The (V,P) value at which the processes intersect,
        #Values for each process interpolated between start and end point
        self.processes = [] 
        self.intersects = []
        self.values = []
        
        #The temperature at each process' start point
        #The work done on the system during each process
        #The heat added to the system during each process
        #and the change in energy during each process
        self.T_values = [] 
        self.W_values = []
        self.Q_values = []
        self.delta_E_values = []

        
        if random_default:
            #generate 4 random PV curves as default
            self.gen_4step_cycle()
            #generate T and n for each step based on ideal gas law
            self.set_4step_T_and_n()
            self.calc_QWE()

    def get_VP_CW(self):
        """returns all (V,P) pairs in clockwise order, 
        starting with the lower right hand point
        """
        return [(self.V0,self.P0), self.intersects[0],
                (self.Vf,self.Pf), self.intersects[1]]

    def calc_QWE(self):
        """Calculates Q, W, and delta E for each step"""
        
        #get V,P and T at each step, repeating the initial value at the end
        #for loop completion
        all_VP = self.get_VP_CW()
        all_VP.append(all_VP[0])
        T_vals = list(self.T_values) #copy array
        T_vals.append(T_vals[0])
        
        #get the process name at each step
        #flatten the process name tuple
        flat_p = [self.processes[0][0],self.processes[0][1],
                  self.processes[1][1],self.processes[1][0]]
        #loop through the (V,P) pairs, T values, and process names, computing
        #the Q, W, and E using the appropriate equations
        for (V0,P0),(Vf,Pf),T0,Tf,pname in zip(all_VP[:-1],all_VP[1:],
                                               T_vals[:-1],T_vals[1:],
                                               flat_p):
            delta_T = Tf - T0
            
            #For all processes: delta_E = nCv*delta_T
            if pname == "T":
                delta_E = 0.
            else:
                delta_E = self.n * intersect.Cv * delta_T            
            
            #Isochoric process: W = 0
            if pname == "C":            
                W = 0.    
            #Isobaric process: W = -PdV = -nRdT
            if pname == "B":
                W = -intersect.R * self.n * delta_T
            #Adiabatic process: W = delta_E                
            if pname == "A":
                W = delta_E
            #Isothermal process: W = -nRTln(Vf/V0)
            if pname == "T": 
                W = -self.n*intersect.R*T0*np.log(Vf/V0)
                
            #For all processes: Q = dE - W
            Q = delta_E - W
            
            self.Q_values.append(Q)
            self.W_values.append(W)
            self.delta_E_values.append(delta_E)        
        
    def gen_4step_cycle(self):
        """Creates a four step complete thermal cycle"""
        start = ['A','B','C','T']
        end = ['A','B','C','T']
        valids = np.zeros(2,dtype = bool)
        
        #generate two valid 2 step processes from (V0,P0) to (Vf,Pf)
        for i in range(2):
            while not valids[i]:
                process_start = start[np.random.randint(len(start))]
                process_end = end[np.random.randint(len(end))]
                try:
                    Vm,Pm = intersect.find_intersect(self.V0,self.P0,
                                                     self.Vf,self.Pf,
                                                     process_start,process_end)
                except TypeError:
                    Vm = None
                    Pm = None
                    
                if Vm is not None and Pm is not None:
                    #print process_start,process_end,Vm,Pm
                    valids[i] = True
                    #write info about the processes
                    if process_start != 'C':
                        eqn_str = intersect.eqn_str(self.V0,self.P0,process_start)
                        p1_V,p1_P = solve_eq.valsfstr(eqn_str, Vm, self.V0)
                    else:
                        p1_V = np.array([Vm,Vm])
                        p1_P = np.array([self.P0,Pm])
                    
                    if process_end != 'C':
                        eqn_str = intersect.eqn_str(self.Vf,self.Pf,process_end)
                        p2_V,p2_P = solve_eq.valsfstr(eqn_str, self.Vf, Vm)
                    else:
                        p2_V = np.array([Vm,Vm])
                        p2_P = np.array([self.Pf,Pm])
                
                    self.values.append((p1_V,p1_P,p2_V,p2_P))
                    self.processes.append((process_start,process_end))
                    self.intersects.append((Vm,Pm))
                    
                    start.remove(process_start)
                    end.remove(process_end)
        #sort the processes so that the "lower" one is listed first
        inter = self.intersects
        if inter[0][0] * inter[0][1] > inter[1][0] * inter[1][1]:
            self.processes = self.processes[::-1]            
            self.intersects = self.intersects[::-1]
            self.values = self.values[::-1]            
        

    def set_4step_T_and_n(self,T0 = None,min_T0 = 100, max_T0 = 1000):
        """sets an arbitrary T value for the initial point 
        and then calculates n based on the ideal gas law
        """
        
        #set initial T value
        if T0 is None:
            self.T_values.append(np.random.randint(min_T0,max_T0))
        else:
            self.T_values.append(T0)
        
        #n = PV/RT
        self.n = self.P0 * self.V0 / (intersect.R * self.T_values[0])

        all_VP = self.get_VP_CW()
        
        #Tf = T0PfVf/(P0V0)
        for (V0,P0),(Vf,Pf) in zip(all_VP[:-1],all_VP[1:]):
            T0 = self.T_values[-1]
            self.T_values.append(T0*Pf*Vf/(P0*V0))

    def calc_efficiency(self):
        #n = 1 - Qc/Qh
        Qh = 0.
        Qc = 0.
        for Q in self.Q_values:
            if Q > 0:
                Qh += Q
            else:
                Qc += -Q
                
        return 1.-Qc/Qh
        
    def get_processes(self):
        return self.processes,self.intersects,self.values

                
    def draw(self, given_value = 'T0'):
        from matplotlib import pyplot as plt
        """Plot points in cycle"""
        #color for each type of process
        process_colors = {"A":"r","B":"g","C":"m","T":"b"}
        p,i,v = self.get_processes()

        #plot the processes        
        plt.plot(v[0][0],v[0][1],color = process_colors[p[0][0]],lw = 3)
        plt.plot(v[0][2],v[0][3],color = process_colors[p[0][1]],lw = 3)
        plt.plot(v[1][0],v[1][1],color = process_colors[p[1][0]],lw = 3)
        plt.plot(v[1][2],v[1][3],color = process_colors[p[1][1]],lw = 3)
        
        #plot the intersection points
        plt.plot(self.V0,self.P0,'c<')
        plt.plot(self.Vf,self.Pf,'c>')
        plt.plot(i[0][0],i[0][1],'c^')
        plt.plot(i[1][0],i[1][1],'cv')

        #set tick labels
        vol_ticks = [self.V0]
        [vol_ticks.append(i[vol][0]) for vol in range(len(i))]
        vol_ticks.append(self.Vf)
        vol_ticks = list(set(vol_ticks))
        
        
        pres_ticks = [self.P0]                
        [pres_ticks.append(i[pres][1]) for pres in range(len(i))]
        pres_ticks.append(self.Pf)
        pres_ticks = list(set(pres_ticks))
        
        vol_lbls = ["%.2f"%vol for vol in vol_ticks]        
        pres_lbls = ["%.2e"%pres for pres in pres_ticks]
        
        plt.xticks(vol_ticks,vol_lbls,rotation=45)
        plt.yticks(pres_ticks,pres_lbls,rotation = 45)

        #default axis labels and title        
        plt.xlabel(r"$\mathdefault{Volume (m^{3})}$")
        plt.ylabel("Pressure (Pa)")
        plt.title("4 Step Heat Engine")
        
        
        dV = (self.V0-self.Vf)/5.
        dP = (self.Pf-self.P0)/5.
        
        plt.gca().set_xlim([self.Vf-dV,self.V0+dV])
        plt.gca().set_ylim([self.P0-dP,self.Pf+dP])        
        #label the given value at the initial point
        if given_value == 'T0':
            plt.text(self.V0-dV/5.,self.P0-dP/2.,
                     r"$\mathdefault{T_0= %sK}$"%self.T_values[0],
                     color = 'r',fontweight='bold')
        elif given_value == 'n':
            plt.text((self.V0+self.Vf)/2.,self.P0-dP/1.5,
                     "n = %.2e mol"%self.n,
                     color = 'b', fontweight = 'bold', ha = 'center')
            
            
def far_apart(V0,P0,Vf,Pf,dV,dP):
    #checks that two points are significantly far from each other
    V_diff = np.abs(Vf-V0) >= dV 
    P_diff = np.abs(P0-Pf) >= dP
    #or that they will be truncated to the same number when displayed
    P_look_same = "%.2e"%P0 == "%.2e"%Pf
    V_look_same = "%.2f"%V0 == "%.2f"%Vf
    return (V_diff or V_look_same) and (P_diff or P_look_same)

def plot_heat_cycle(HC, fname = None, figsize=(8,6),given_value = 'T0'):
    from matplotlib import pyplot as plt    
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot(111)
    HC.draw(given_value = given_value)
    ax.grid()
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname,bbox_inches = "tight")

def text_QWE_table(HC,fname):
    T_vals = list(HC.T_values)
    T_vals.append(T_vals[0])
    with open(fname,'w') as f:
        for i,Q,W,E,T0,Tf in zip(range(len(HC.Q_values)),HC.Q_values,
                         HC.W_values,HC.delta_E_values,
                         T_vals[:-1],T_vals[1:]):
            j = i+2
            if j > len(HC.Q_values): j = 1
            f.write("%s%s %.2e %.2e %.2e %.2e %.2e\n"%(i+1,j,Q,W,E,T0,Tf))
def gen_img_and_png(img_fname = None, table_fname = 'default.txt'):
    #from matplotlib import pyplot as plt
    """test case"""
    #makes sure the midpoints for the processes are not too close to the end
    #points    
    readable_graph = False
    while not readable_graph:
        readable_graph = True

        V0 = np.random.randint(1,40)/4.
        Pf = np.random.rand()*10**np.random.randint(1,5)
        P0 = Pf * .75*np.random.rand()
        Vf = V0 * .75*np.random.rand()
        
        dV = (V0-Vf)/10.
        dP = (Pf-P0)/10.        
        
        
        hc = HeatCycle(V0,P0,Vf,Pf)
        #hc.gen_4step_cycle()
        p,i,v = hc.get_processes()
        points =  ((V0,P0),(Vf,Pf),i[0],i[1])
        
        #make sure none of the points are too close to eachother
        for p1 in points:
            for p2 in points:
                if not far_apart(p1[0],p1[1],p2[0],p2[1],dV,dP):
                    readable_graph = False
                    
                        
    plot_heat_cycle(hc,fname = img_fname, given_value=['n','T0'].pop(np.random.randint(2)))
    text_QWE_table(hc,table_fname)
    
    
if __name__ == '__main__':
    import sys
    #python heat_cycle.py image_fname table_fname
    gen_img_and_png(sys.argv[1],sys.argv[2])
                   