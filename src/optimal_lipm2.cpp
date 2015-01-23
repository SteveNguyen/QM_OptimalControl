
#include "Models/Lipm2/Lipm2Description.h"


#include <time.h>
using namespace std;



int main(int argc, char*argv[])
{




   time_t start, stop;
   double diff;
   Lipm2Manager F;
    //LipmOneGaussManager F;
    //LipmLearningManager F;


       //F.build(DRAWN_POLICY, INIT_BETA);

/** Objective specification (now specified in *_config.h files in include/Models) :  ***/

//    vector<double> obj;
//    obj.push_back(0.0); // position
//    obj.push_back(0.5); // speed
//    F.setObjective(obj);


/** Optimization **/

    if (argc>1 && !strcmp(argv[1],"-l"))
    {
        F.build(DRAWN_POLICY, INIT_BETA);
        F.load("BE.ser");
    }
    else
    {
        start = time(NULL);
        F.build(DRAWN_POLICY, INIT_BETA);
        F.compute();
        stop = time(NULL);
        diff = difftime(stop, start);
        printf("Elapsed time %f s\n", diff);

        
        F.save("BE.ser");
    }



    F.toPlot();


/*** Control simulation: ***/

        // F.run();

    return 0;
}
