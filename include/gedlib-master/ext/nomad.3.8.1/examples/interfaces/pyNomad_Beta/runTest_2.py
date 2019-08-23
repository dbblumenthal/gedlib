import PyNomad
import sys

def bb(x):
    try:
        dim = x.get_n()
        f = x.get_coord(4)
        x.set_bb_output(0,f)
        g1 = sum([(x.get_coord(i)-1)**2 for i in xrange(dim)])-25
        x.set_bb_output(1,g1)
        g2 = 25-sum([(x.get_coord(i)+1)**2 for i in xrange(dim)])
        x.set_bb_output(2,g2)
        # test = 1* (1/0) # ---> uncomment to trigger an exception
    except:
        print "Unexpected error:", sys.exc_info()[0]
        return -1
    return 1

ub = [ 5,6,7,10,10]

params = ['BB_OUTPUT_TYPE OBJ PB EB','MAX_BB_EVAL 100','X0 * 0' , 'LOWER_BOUND * -6' ] 
[ x_return , f_return , h_return, nb_evals , nb_iters ,  stopflag ] = PyNomad.optimize(bb,[],[],ub,params)
print ('\n NOMAD outputs \n X_sol={} \n F_sol={} \n H_sol={} \n NB_evals={} \n NB_iters={} \n'.format(x_return,f_return,h_return,nb_evals,nb_iters))
