#ifndef CONTROLLER_H
#define CONTROLLER_H

//#include "Simulator.h"
//#include "Optimizer.h"
#include "Learner.h"

using namespace std;

class Controller
{
    public:
        Controller(Optimizer& eng_, ModelLearner& learner_);

        virtual Value& control(Value state_val, Value& obj)
        {

            return eng.chooseAction(state_val, obj);
        }
        virtual void run() {;}
        virtual ~Controller();
    protected:
        Optimizer& eng;
        ModelLearner& learner;
        Variable& state;
        Variable& action;
    private:
};

class LearningController : public Controller
{
    public:
        LearningController(Optimizer& eng_, ModelLearner& learner_, Simulator& sim_ /*actually irrelevant*/) : Controller(eng_, learner_), learner(action, 0, state /*actually empty*/, sim_, 0), u_d(action), x_d(state), u_c(action), x_c(state)
        {
            x_d.reset();
            x_c.continuize(x_d);   // default right value for the learner
            Distributions::uniform(hist_U, action, 0);
            gen.seed(time(NULL));
        }
        Value& control(Value state_val, Value& obj)
        {

                //with opposite hist
            Distributions::Distrib1D opp_hist_U(Distributions::opposite(hist_U, action, 0));
            u_d[0]=Distributions::draw(Distributions::toVect(opp_hist_U, action, 0), gen);




            /* //random */
            /* Distributions::Distrib1D uni; */
            /* Distributions::uniform(uni, action, 0); */
            /* u_d[0]=Distributions::draw(Distributions::toVect(uni, action, 0), gen); */
            
            u_c.continuize(u_d);
            learner.addPoint(u_c,x_c);
            hist_U=learner.getDistribution(x_d);
            return u_d;
        }
    private:
        CndHistogram1DLearner learner;
        Distributions::Distrib1D hist_U;
        Value u_d, x_d;
        ContValue u_c, x_c;
        boost::mt19937 gen;
        int u;

};

#endif // CONTROLLER_H
