#ifndef LEG2DHDESCRIPTION_H
#define LEG2DHDESCRIPTION_H

#include "Specification.h"
#include "Manager.h"


#include "Config.h"



class Leg2DHDescription : public ModelDescription
{
    public:
        Leg2DHDescription(Variable& state, Variable& action, Simulator& sim) : ModelDescription(state,action,sim)
        {
            ;
        }

    #ifndef _COST_ACC
    double cost(const Value& state_, const Value& action_)
        {
            return 1.0;  //abs(action_[0]-(action.cardinality()-1)/2.0)+1;
                
        }
    #else
    
    double cost(const Value& state_, const Value& action_)
    {
        Value s(state_);
        Value a(action_);
        
        ContValue X_cval(state);
        ContValue U_cval(action);       
        
        X_cval.continuize(s);
        U_cval.continuize(a);

        double cost=1.0;

        cost+=X_cval[2]+pow(U_cval[0],2)+pow(U_cval[1],2);
        
        return cost;
        
        
    }
    #endif


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

        ~Leg2DHDescription();
    protected:

};


class Leg2DHSimulator : public Simulator
{
    public:
       // using
        ContValue& simulate(double dt)
        {
            return Simulator::simulate(dt);

        }



        void simulate(ContValue& state, ContValue& action, double dt);
            /*
        void setFriction(double f_) {f=f_;}
        void setSlope(double s) {slope=s;}
        Leg2DHSimulator(Variable& state, Variable& action, double defaultDT_, double min_=0.0, double max_=2*M_PI) : Simulator(state,action,defaultDT_), min(min_), max(max_), f(0.0), slope(0.0) {;}
            */
  Leg2DHSimulator(Variable& state, Variable& action, double defaultDT_) : Simulator(state,action,defaultDT_) {;}
    protected:
//        double min, max, f, slope;


};

class Leg2DHOneGaussDescription : public Leg2DHDescription
{
    public:
        Leg2DHOneGaussDescription(Variable& state, Variable& action, Simulator& sim) : Leg2DHDescription(state, action, sim) {;}
        void getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U);
};

//class Leg2DHOneGaussSimulator : public Leg2DHSimulator
//{
//    Leg2DHOneGaussSimulator(Variable& state, Variable& action, double min_=0.0, double max_=2*M_PI) : Leg2DHSimulator(state, action, min_, max_) {;}
//    void simulate(ContValue& state, ContValue& action, double dt);
//};

class Leg2DHManager : public Manager
{
    public:
       // Leg2DHManager();
        void buildVariables();
        void buildSimulator();
        void buildModel();
        void buildOptimizer(bool drawnPolicy, double initBeta);
        void buildDistribLearner();
        void buildModelLearner();
        void buildController();
        void buildObjective();
        void run();
        virtual ~Leg2DHManager();
    protected:
    private:

};

class Leg2DHOneGaussManager : public Leg2DHManager
{
    public:
       // Leg2DHManager();
        void buildModel()
        {
            model = new Leg2DHOneGaussDescription(*state, *action, *sim);
        }

        void buildDistribLearner()
        {
            distrib1DLearners.push_back(new CndHistogram1DLearner(*state, 0, (*state), *sim, 0));
            //distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 0, (*state)^(*action), NB_SIG, *sim, SIGMA_P, 4));
            distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 1, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
        }

        void buildModelLearner()
        {
            learner = new ModelLearnerOneGaussOneDirac(*optimizer, distrib1DLearners, *state, (*state)^(*action), *sim, DT);
        }
//        void buildDistribLearner();
//        void buildModelLearner();

//        virtual ~Leg2DHOneGaussManager() {;}
    protected:
    private:

};

class Leg2DHLearningManager : public Leg2DHManager
{
    public:
       // Leg2DHManager();
        void buildController()
        {
            controller = new LearningController(*optimizer, *learner, *sim);
        }


//        void buildDistribLearner();
//        void buildModelLearner();

//        virtual ~Leg2DHOneGaussManager() {;}
    protected:
    private:

};


#endif // PENDULUMDESCRIPTION_H
