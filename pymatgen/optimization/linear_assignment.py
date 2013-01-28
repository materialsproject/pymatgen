import numpy as np

"""
Implementation of the LAPJV algorithm described in:

R. Jonker, A. Volgenant. A Shortest Augmenting Path Algorithm for 
Dense and Sparse Linear Assignment Problems. Computing 38, 325-340
(1987)

Uses numpy where possible
Minimizes the cost of the assignment
"""

class LinearAssignment():
    
    def __init__(self, costs):
        self.c = np.array(costs)
        self.n = len(costs)
        
        self._x = np.zeros(self.n, dtype = np.int)-1
        self._y = np.zeros(self.n, dtype = np.int)-1
        
        #preprocess
        if self._column_reduction():
            #only augment if column reduction doesn't find full soln
            #dual row and column variables
            self._augmenting_row_reduction()
            #initialize the reduced costs
            self._update_cred()
            while np.min(self._x) < 0:
                self._augment()
            
        self.solution = self._x
        self._min_cost = None
    
    
    @property
    def min_cost(self):
        if self._min_cost:
            return self._min_cost
        
        self._min_cost = np.sum(self.c[np.arange(len(self.c)), self.solution])
        return self._min_cost
    
    
    def _column_reduction(self):
        #column_reduction
        i1, j = np.unique(np.argmin(self.c, axis = 0), return_index = True)
        self._x[i1] = j
        
        if len(i1) == self.n:
            #problem is solved, return
            return False
        
        self._y[j] = i1
        
        #reduction_transfer
        #tempc is array with previously assigned entry masked
        self._v = np.min(self.c, axis = 0)
        tempc = self.c.copy()
        tempc[i1, j] = np.max(tempc.flatten()) * 10
        mu = np.min(tempc[i1, :] - self._v[None, :], axis = 1)
        self._v[j] -= mu
        return True
        
        
    def _augmenting_row_reduction(self):
        unassigned = np.where(self._x == -1)[0]
        for i in unassigned:
            while True:
                #find smallest 2 values and indices
                temp = self.c[i] - self._v
                j1 = np.argmin(temp)
                u1 = temp[j1]
                temp[j1] = np.max(temp) + 1
                j2 = np.argmin(temp)
                u2 = temp[j2]
                
                if u1 < u2:
                    self._v[j1] -= u2-u1
                elif self._y[j1] != -1:
                    j1 = j2
                k = self._y[j1] 
                if k != -1:
                    self._x[k] = -1
                    self._x[i] = j1
                    self._y[j1] = i
                    i = k
                if u1 == u2 or k == -1:
                    break
                    
    
    def _update_cred(self):
        """
        updates the reduced costs
        """
        ui = np.diag(self.c[:, self._x]) - self._v[self._x]
        self.cred = self.c - ui[:, None] - self._v[None, :]
        
    
    def _augment(self):
        """
        Finds a minimum path and adds it to the matching
        """
        self._build_tree()
        #update prices
        delta = self._d[self._ready] - self._mu
        self._v[self._ready] += delta
        
        #augmentation
        while True:
            self._i = self._pred[self._j]
            self._y[self._j] = self._i
            k = self._j
            self._j = self._x[self._i]
            self._x[self._i] = k
            if self._i == self._istar:
                break
        self._update_cred()
    
    
    def _build_tree(self):
        """
        builds the tree finding an augmenting path
        """
        #find unassigned i*
        self._istar = np.argmin(self._x)
        
        #compute distances
        self._d = self.c[self._istar] - self._v
        self._pred = np.zeros(self.n, dtype = np.int) + self._istar
        
        self._ready = np.zeros(self.n, dtype = np.bool)
        self._scan = np.zeros(self.n, dtype = np.bool)
        self._todo = np.zeros(self.n, dtype = np.bool) + True

        while True:
            #populate scan with minimum reduced distances
            if np.max(self._scan) == 0:
                self._mu = np.min(self._d[self._todo])
                self._scan[np.where(self._d == self._mu)] = 1
                self._todo[self._scan] = 0
                if np.min(self._y * self._scan) < 0:
                    self._j = np.argmin(self._y * self._scan)
                    return
                
            #pick jstar from scan (scan always has at least 1)
            self._jstar = np.argmax(self._scan)
            
            #pick i associated with jstar
            self._i = self._y[self._jstar]
            
            self._scan[self._jstar] = 0
            self._ready[self._jstar] = 1
            
            #find shorter distances
            newdists = self._mu + self.cred[self._i, :]
            shorter = (newdists < self._d) * self._todo
            
            #update distances
            self._d[shorter] = newdists[shorter]
            
            #update pred
            self._pred[shorter] = self._i
            
            for self._j in np.argwhere((self._d == self._mu) * self._todo).flatten():
                if self._y[self._j] == -1:
                    return
                self._scan[self._j] = 1
                self._todo[self._j] = 0
