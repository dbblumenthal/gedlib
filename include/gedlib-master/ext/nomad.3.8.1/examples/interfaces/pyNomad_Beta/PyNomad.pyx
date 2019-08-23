# C.Tribes sept. 6th, 2016 --- PyNomad: Beta integration of Nomad in Python 
# PyNomad 1.0

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp cimport bool

from cython.operator cimport dereference as deref, preincrement as inc
from multiprocessing import Process, Queue

def version():
    printPyNomadVersion()
 
# Define the interface function to display nomad general information
def usage():
    printPyNomadUsage()
 
# Define the interface function to display nomad general information
def info():
    printPyNomadUsage()
    printPyNomadInfo()

# Define the interface function to get nomad help
def help(about=''):
    printNomadHelp(about)
    
def __doc__():
    cdef string about;
    printPyNomadUsage()
    help(about)	  

# Define the interface function to perform optimization          
def optimize(f , pX0, pLB,pUB ,params):
    cdef PyNomadEval_Point u_feas = PyNomadEval_Point()
    cdef PyNomadEval_Point u_infeas = PyNomadEval_Point()
    cdef int run_status
    cdef int nb_evals = 0
    cdef int nb_iters = 0
    cdef double f_return
    cdef double h_return
    x_return=[]
    run_status = runNomad(cb, cbL, <void*> f, <vector[double]&> pX0, <vector[double]&> pLB, <vector[double]&> pUB , <vector[string]&> params , u_feas.c_ep ,u_infeas.c_ep , nb_evals, nb_iters)
    if u_feas.c_ep != NULL:
         f_return = u_feas.get_f()
         h_return = 0
         for i in xrange(u_feas.get_n()):
             x_return.append(u_feas.get_coord(i))
         del u_feas.c_ep  
    if u_infeas.c_ep != NULL:
         f_return = u_infeas.get_f()
         h_return = u_infeas.get_h()
         for i in xrange(u_infeas.get_n()):
             x_return.append(u_infeas.get_coord(i))
         del u_infeas.c_ep      
    return [ x_return , f_return, h_return, nb_evals , nb_iters  , run_status ]
   
cdef extern from "Double.hpp" namespace "NOMAD":
    cdef cppclass Double:
        const double & value()
        bool is_defined()
    
cdef class PyNomadDouble:
    cdef Double c_d 
    def value(self):
        return self.c_d.value()
    def is_defined(self):
        return self.c_d.is_defined()
	       
cdef extern from "Eval_Point.hpp" namespace "NOMAD":
    cdef cppclass Eval_Point:
        const double & value(int i)      
        const Double & get_f()
        const Double & get_h()
        void set_bb_output(int i, double & v)
        int get_n()
        int get_m()
                     
cdef class PyNomadEval_Point:
    cdef Eval_Point *c_ep 
    def get_coord(self, int i):       
        return self.c_ep.value(i)
    def set_bb_output(self,int i, double v):
        self.c_ep.set_bb_output(i,v)
    def get_f(self):
        cdef PyNomadDouble f = PyNomadDouble()
        f.c_d=self.c_ep.get_f()
        cdef double f_d
        if ( f.is_defined() ): 
            f_d = f.value()
        else:
            f_d = float('inf')
        return f_d  
    def get_h(self):
        cdef PyNomadDouble h = PyNomadDouble()
        h.c_d=self.c_ep.get_h()
        cdef double h_d
        if ( h.is_defined() ):
            h_d = h.value()
        else:
            h_d = 0
        return h_d  

    def get_n(self):
        cdef int n
        n = self.c_ep.get_n()
        return n
    def get_m(self):
        cdef int m
        m = self.c_ep.get_m()
        return m

cdef extern from "nomadCySimpleInterface.cpp":
    ctypedef int (*Callback)(void * apply, Eval_Point& x)
    ctypedef int (*CallbackL)(void * apply, list[Eval_Point *] & x)
    void printPyNomadInfo()
    void printPyNomadUsage()
    void printNomadHelp( string about)
    void printPyNomadVersion()
    int runNomad(Callback cb, CallbackL cbL, void* apply, vector[double] &X0, vector[double] &LB, vector[double] &UB , vector[string] & params , Eval_Point *& best_feas_sol ,Eval_Point *& best_infeas_sol , int & nb_evals, int & nb_iters) except+

# Define callback function for a single Eval_Point ---> link with Python     
cdef int cb(void *f, Eval_Point & x):
     cdef PyNomadEval_Point u = PyNomadEval_Point()
     u.c_ep = &x
     return (<object>f)(u)  
 
# Define callback function for block evaluation of a list of Eval_Points ---> link with Python    
cdef int cbL(void *f, list[Eval_Point *] & x):
      cdef size_t size = x.size()
      cdef PyNomadEval_Point u
      cdef list[Eval_Point *].iterator it = x.begin()
      cdef Eval_Point *c_ep
      
      # Start the process for each Eval_Point of the list
      out = Queue()
      for i in xrange(size):
         u = PyNomadEval_Point()
         c_ep = deref(it)
         u.c_ep = c_ep
         proc = []
         p = Process(target=<object>f, args=(u,out,)) # u is copied
         p.start()
         proc.append(p)
         inc(it)
       
      # Wait for all the processes to end         
      for p in proc:
         p.join()
        
      # Update the Eval_Points bb_output from the out Queue  
      it = x.begin()
      for i in xrange(size):
          c_ep = deref(it)
          bb_out = out.get()
          if type(bb_out) != type(float) or type(bb_out) != type(int):
              for j in xrange(len(bb_out)):
                  c_ep.set_bb_output(j,bb_out[j])
          else:
              c_ep.set_bb_output(0,bb_out)
                          
          inc(it)  
   
      return 1   # 1 is success
 
