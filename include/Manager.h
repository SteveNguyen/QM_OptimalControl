#ifndef MANAGER_H
#define MANAGER_H

#include "Controller.h"

class Manager
{
    public:
  //      SimulatorManager();
        virtual void buildVariables()=0;
        virtual void buildSimulator()=0;
        virtual void buildModel()=0;
        virtual void buildOptimizer(bool drawnPolicy, double initBeta)=0;
        virtual void buildDistribLearner()=0;
        virtual void buildModelLearner()=0;
        virtual void buildController()=0;
        virtual void buildObjective()=0;
        void build(bool drawnPolicy, double initBeta)
        {
            buildVariables();
            buildSimulator();
            buildModel();
            buildOptimizer(drawnPolicy, initBeta);
            buildDistribLearner();
            buildModelLearner();
            buildController();
            obj = new Value(*state);
            obj_c = new ContValue(*state);
            buildObjective();
            obj->discretize(*obj_c);

        }
        virtual void run()=0;
        Simulator& getSimulator() {return *sim;}
        Controller& getController() {return *controller;}
        ModelLearner& getLearner() {return *learner;}
        Optimizer& getOptimizer() {return *optimizer;}
        void setObjective(vector<double>& objVect)
        {
            for (unsigned int i=0; i<state->nbDim(); i++)
                (*obj_c)[i]=objVect[i];
            obj->discretize(*obj_c);
            cout<<"Objective = (";
            for (unsigned int i=0; i<state->nbDim(); i++)
                cout<<(*obj)[i]<<" ";
            cout<<")"<<endl;
        }
        void load(string file)
        {
            optimizer->load(file);
        }
        void save(string file)
        {
            optimizer->save(file);
        }
        void compute()
        {
            optimizer->computeLandscape(false, *obj);
            optimizer->computePolicy(*obj);
        }

        void computeT() //compute for the transposed graph (single source problem)
        //not very efficient since it first transpose the multi source graph
        {
            optimizer->computeLandscapeT(false, *obj);
            optimizer->computePolicy(*obj);
        }
        
        void toPlot()
        {
            optimizer->toPlot(*obj);
        }
        void setDir(string dir)
        {
            optimizer->setDir(dir);
        }
  //      virtual ~SimulatorManager();
    protected:
        Variable* state;
        Variable* action;
        Value* obj;
        ContValue* obj_c;
        Simulator* sim;
        ModelDescription* model;
        Optimizer* optimizer;
        vector<CndDistribution1DLearner*> distrib1DLearners;
        ModelLearner* learner;
        Controller* controller;


    private:


};

#endif // MANAGER_H
