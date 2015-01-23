#ifndef PENDULUMDESCRIPTION_H
#define PENDULUMDESCRIPTION_H

#include "Specification.h"
#include "Manager.h"


#include "TodorovConfig.h"
//#include "WalkingConfig.h"
//#include "../../../../optimal_walker/include/RobotConfig.h"


class PendulumDescription : public ModelDescription
{
    public:
        PendulumDescription(Variable& state, Variable& action, Simulator& sim) : ModelDescription(state,action,sim)
        {
            ;
        }

        
        double cost(const Value& state_, const Value& action_)
        {
                //return 1.0;  //abs(action_[0]-(action.cardinality()-1)/2.0)+1;

            if(state_[0]==XCARD/2 && state_[1]==VCARD/2 && action_[0]==UCARD/2)
                return 0.0;

            else
                return 1.0;
            
            
        }
        


        /*
    double cost(const Value& state_, const Value& action_)
        {
                //return 1.0;  //abs(action_[0]-(action.cardinality()-1)/2.0)+1;


            Value s(state_);
            Value a(action_);
            
            ContValue X_cval(state);
            ContValue U_cval(action);
            
       
            
            X_cval.continuize(s);
            U_cval.continuize(a);

            
            
            if(state_[0]==XCARD/2 && state_[1]==VCARD/2 && action_[0]==UCARD/2)
                return 0.0;

            else
                return (fabs(X_cval[0]))+1.0;//+fabs(U_cval[0]));
                    //return 1.0/(1.0-3.9*fabs(U_cval[0]));
                    //return (fabs(X_cval[0])+100*fabs(U_cval[0]));
                
            
        }
        */
        void getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U);

        //vector<double> learnParams(Value& Xt, Value& U, unsigned int nb_sample_per_dim);

        ~PendulumDescription();
    protected:

};


class PendulumSimulator : public Simulator
{
    public:
       // using
        ContValue& simulate(double dt)
        {
            return Simulator::simulate(dt);

        }



        void simulate(ContValue& state, ContValue& action, double dt);
        void setFriction(double f_) {f=f_;}
        void setSlope(double s) {slope=s;}
        PendulumSimulator(Variable& state, Variable& action, double defaultDT_, double min_=0.0, double max_=2*M_PI) : Simulator(state,action,defaultDT_), min(min_), max(max_), f(0.0), slope(0.0) {;}

    protected:
        double min, max, f, slope;


};

class PendulumOneGaussDescription : public PendulumDescription
{
    public:
        PendulumOneGaussDescription(Variable& state, Variable& action, Simulator& sim) : PendulumDescription(state, action, sim) {;}
        void getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U);
};

//class PendulumOneGaussSimulator : public PendulumSimulator
//{
//    PendulumOneGaussSimulator(Variable& state, Variable& action, double min_=0.0, double max_=2*M_PI) : PendulumSimulator(state, action, min_, max_) {;}
//    void simulate(ContValue& state, ContValue& action, double dt);
//};

class PendulumManager : public Manager
{
    public:
       // PendulumManager();
        void buildVariables();
        void buildSimulator();
        void buildModel();
        void buildOptimizer(bool drawnPolicy, double initBeta);
        void buildDistribLearner();
        void buildModelLearner();
        void buildController();
        void buildObjective();
        void run();
        virtual ~PendulumManager();
    protected:
    private:

};

class PendulumOneGaussManager : public PendulumManager
{
    public:
       // PendulumManager();
        void buildModel()
        {
            model = new PendulumOneGaussDescription(*state, *action, *sim);
        }

        void buildDistribLearner()
        {
            distrib1DLearners.push_back(new CndHistogram1DLearner(*state, 0, (*state), *sim, 0));
            
            //distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 0, (*state)^(*action), NB_SIG, *sim, SIGMA_P, 4));
            /* distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 0, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT)); */
            distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 1, (*state)^(*action), NB_SIG, *sim, SIGMA_S, LEARN_INIT_WEIGHT));
        }

        void buildModelLearner()
        {
            learner = new ModelLearnerOneGaussOneDirac(*optimizer, distrib1DLearners, *state, (*state)^(*action), *sim, DT);
        }
//        void buildDistribLearner();
//        void buildModelLearner();

//        virtual ~PendulumOneGaussManager() {;}
    protected:
    private:

};

class PendulumLearningManager : public PendulumManager
{
    public:
       // PendulumManager();
        void buildController()
        {
            controller = new LearningController(*optimizer, *learner, *sim);
        }


//        void buildDistribLearner();
//        void buildModelLearner();

//        virtual ~PendulumOneGaussManager() {;}
    protected:
    private:

};


#endif // PENDULUMDESCRIPTION_H
