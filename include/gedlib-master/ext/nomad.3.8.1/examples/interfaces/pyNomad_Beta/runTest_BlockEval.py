import PyNomad
import sys

# This example of blackbox function is for multiprocess: several bb are called in parallel
# Multiprocess blackbox evaluations request to provide the BB_MAX_BLOCK_SIZE in NOMAD params
# 
# In this case blackbox function requires 2 arguments
# The first argument x is similar to a NOMAD::Eval_Point --> access the coordinates with get_coord() function
# The blackbox output must be put in the second argument bb_out and the bb outputs must be put in the same order as defined in params
# bb_out is a queue to work in a multiprocess.
def bb(x,bb_out):
    try:
        dim = x.get_n()
        f = x.get_coord(4)
        g1 = sum([(x.get_coord(i)-1)**2 for i in xrange(dim)])-25
        g2 = 25-sum([(x.get_coord(i)+1)**2 for i in xrange(dim)])
        bb_out.put([f,g1,g2])
    except:
        print "Unexpected error:", sys.exc_info()[0]
	sys.exit(0)        

x0 = [0,0,0.71,0.51,0.51]
lb = [-6,-6,-6,-6,-6]
ub = [ 5,6,7,10,10]

params = ['BB_OUTPUT_TYPE OBJ PB EB','MAX_BB_EVAL 100' , 'BB_MAX_BLOCK_SIZE 8','DISPLAY_STATS BBE BLK_EVA OBJ'] 
[ x_return , f_return , h_return, nb_evals , nb_iters ,  stopflag ] = PyNomad.optimize(bb,x0,lb,ub,params)
print ('\n NOMAD outputs \n X_sol={} \n F_sol={} \n H_sol={} \n NB_evals={} \n NB_iters={} \n'.format(x_return,f_return,h_return,nb_evals,nb_iters))
