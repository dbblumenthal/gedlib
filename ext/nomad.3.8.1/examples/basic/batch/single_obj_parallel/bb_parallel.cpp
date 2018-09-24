#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <pthread.h>
#include <semaphore.h>

using namespace std;

// Number of threads to be used in parallel
#define NUM_THREADS    4

// A semaphore to manage the threads
sem_t mySemaphore;

static void * eval_x ( void * dummy_eval_arg   )
{
    
    double c1 = 0.0 , c2 = 0.0 ,f;
    double * x = static_cast<double *>(dummy_eval_arg);
    for ( int i = 0 ; i < 5 ; i++)
    {
        c1 += pow ( x[i]-1 , 2 );
        c2 += pow ( x[i]+1 , 2 );
    }
    x[5]= x[4];
    x[6] = c1 - 25;
    x[7] = 25 - c2;
    
    // cout << x[0] << " " << x[1] <<" " <<  x[2] <<" "<< x[3] <<" "<< x[4] << " " << x[5] << " " << x[6] << " " << x[7] << endl;
    
    pthread_exit(NULL);
    
    // Increase the semaphore --->  a new thread can be launched
    sem_post(&mySemaphore);
    
}



int main ( int argc , char ** argv )
{
    
    double f = 1e20, c1 = 1e20 , c2 = 1e20;
    std::vector< double * > X;
    
    if ( argc >= 2 )
    {
        
        ifstream in ( argv[1] );
        
        while( ! in.eof() )
        {
            double *x = new double [8];
            for ( int i = 0 ; i < 5 ; i++ )
            in >> x[i];
            
            X.push_back(x);
            
            //char c;
            // in >> c;
            
            // cout << x[0] << " " << x[1] <<endl;
            
        }
        in.close();
        X.pop_back();
        
        int nb_pts=X.size();
        
        // All threads are created
        pthread_t threads[nb_pts];
        
        // Initialize the semaphore with the number of threads available
        sem_init(&mySemaphore,0,NUM_THREADS);
        
        // Another way to use semaphore
        //  sem_t * mySemaphore;
        //  mySemaphore=sem_open("sem",0,NUM_THREADS);
        
        
        // The list of points is evaluated under the control of a semaphore
        int i=0;
        for (std::vector< double * >::iterator it_x=X.begin(); it_x!=X.end(); ++it_x,++i)
        {
            //		cout << (*it_x)[0] << " " << (*it_x)[1] << " " << (*it_x)[2] << endl;
            int rc=pthread_create(&threads[i], NULL, eval_x,(*it_x));
            if (rc)
            {
                cout << "Error:unable to create thread," << rc << endl;
                return false;
            }
            // wait until value of semaphore is greater than 0
            // and decrement the value of semaphore by 1
            sem_wait(&mySemaphore);
            
        }
        
        int ret;
        
        // cout << nb_pts << endl;
        for (i=0; i<nb_pts; ++i)
        {
            // Wait for all the threads to finish
            ret=pthread_join(threads[i],0);
            if (ret!=0)
            {
                perror("pthread join has failed");
                return false;
            }
            cout << (X[i])[5]  << " " << (X[i])[6] << " " << (X[i])[7] << endl;
            delete [](X[i]);
        }
        // Unlink semaphore if the "other" way is used
        // sem_unlink("sem");
        
        sem_destroy(&mySemaphore);
        
    } 
    
    return 0;
    
}	
