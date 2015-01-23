#ifndef SPECIFICATION_H
#define SPECIFICATION_H


#include <vector>
#include <map>
#include <float.h>
#include <math.h>
#include <iostream>
#include <numeric>
#include <fstream>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//#include "ModelDescription.h"



#define EPSI 0.001


using namespace std;

class Variable
{
    public:

        Variable(vector<double>& mins_, vector<double>& maxs_, vector<unsigned int>& cardPerDim_, vector<bool>& isCircular_);
        Variable(const Variable& x ) : nbDim_(x.nbDim_), mins(x.mins), maxs(x.maxs), cardPerDim(x.cardPerDim), isCirc(x.isCirc), card(x.card) {;}
        double getMins(unsigned int dim=0) const {return mins[dim];}
        double getMaxs(unsigned int dim=0) const {return maxs[dim];}
        const vector<unsigned int>& getCardPerDim() const {return cardPerDim;}
        const unsigned int& cardinality() const {return card;}
        const unsigned int& cardinality(unsigned int dim) const {return cardPerDim[dim];}
       
        unsigned int& nbDim() {return nbDim_;}
        const bool isCircular(unsigned int dim) const {return isCirc[dim];}
//        double getMin(unsigned int dim) {return mins[dim];}
//        double getMax(unsigned int dim) {return maxs[dim];}
        Variable operator^(const Variable &other) const;

    private:
        friend class Value;
        friend class ContValue;
        unsigned int nbDim_;
        vector<double> mins;
        vector<double> maxs;
        vector<unsigned int> cardPerDim;
        vector<bool> isCirc;
        unsigned int card;
};

class ContValue; // forward declaration needed for use in Value::discretize(ContValue&) below

class Value
{
    public:
        Value() {/* cout<<"Variable default constructor, probably used in boost graph ..."<<endl */ ;}
        Value(const Value& x) : var(x.var), value(x.value) {;}
        Value(Variable& v): var(&v), value(v.nbDim()) {;}
        Value(vector<unsigned int>& val): value(val) {;}
        void setVariable(Variable& v) {var=&v;value=vector<unsigned int>(v.nbDim());}
        const vector<unsigned int>& getValueVect() const {return value;}
        vector<unsigned int>& getValueVect() {return value;}        // non-const version (able to modify the value vector
        void reset() {value=vector<unsigned int>(var->nbDim_,0);}
        bool next();
        void discretize(ContValue& val);
        unsigned int& operator[](unsigned int ind) {return value[ind];}
        unsigned int operator[](unsigned int ind) const {return value[ind];}
        bool operator<(const Value &other) const
        {
            return (value<other.value);
        }
        Value operator^(const Value &other) const;
        const bool isCircular(unsigned int dim) const 
        {
            return var->isCircular(dim);
        }
        
    private:
        Variable* var;
        vector<unsigned int> value;
        
        friend class boost::serialization::access;
//#ifndef _NOSERIAL
        template <class Archive> void serialize(Archive &ar, const unsigned int version)
        {
            ar & value; // & policy;
        }
//#endif
};

class ContValue
{
    public:
        ContValue(Variable& var_):var(&var_), cvalue(var_.nbDim()) {;}
        ContValue(const ContValue& x) : var(x.var), cvalue(x.cvalue) {;}
        ContValue(Variable& var_, Value val ):var(&var_), cvalue(var_.nbDim())
        {
            continuize(val);
        }

        ContValue(vector<double>& val): cvalue(val) {;}

        void setVariable(Variable& v) {var=&v;cvalue=vector<double>(v.nbDim());}

        ContValue& operator=(const ContValue& other)
        {
            cvalue=other.cvalue;
            return *this;
        }

        double& operator[](unsigned int ind)
        {
            circularize(ind); // maybe not optimal but quite safe: a ContValue is circularized each time it is accessed (and modified).
            return cvalue[ind];
        }

        ContValue operator^(const ContValue &other) const;

        void continuize(Value& val);
        const bool isCircular(unsigned int dim) const 
        {
            return var->isCircular(dim);
        }
        
        void circularize(unsigned int dim);
        vector<double>& getValueVect() {return cvalue;}
        Variable& getVar() {return *var;}
    private:
        Variable* var;
        vector<double> cvalue;
};


class Distributions
{
    public:
        typedef std::map<unsigned int, double> Distrib1D;
        typedef std::map<Value, double> Distrib;
        static void gaussian(Distrib1D& P, Variable& X, unsigned int dim, double mu, double sigma, double nbSig_);
        static void vonMises(Distrib1D& P, Variable& X, unsigned int dim, double mu, double kappa);
        static void dirac(Distrib1D& P, Variable& X, unsigned int dim, double val);
        static void dirac(Distrib& P, Variable& X, unsigned int dim, double val);
        static void uniform(Distrib1D& P, Variable& X, unsigned int dim);

        static bool normalize (vector<double>& v)
        {
            double sum=0.0;
            bool b=true;
            for (vector<double>::iterator it=v.begin(); it<v.end(); it++)
                sum+=*it;
            if (abs(sum)<FLT_EPSILON)
            {
                cerr<<"WARNING: normalizing a null vector returns a uniform distribution"<<endl;
                b=false;
            }

            for (vector<double>::iterator it=v.begin(); it<v.end(); it++)
                if (abs(sum)<FLT_EPSILON)
                {
                    *it=1.0/double(v.size());
                } else
                    *it/=sum;
            return b;
        }
        static vector<double> toVect(Distrib1D& P, Variable& V, unsigned int dim)
        {
            vector<double> res;
           //= Value val(V);
          //  val.reset();
            Distrib1D::iterator it;
            //cout<<"(";
            for (unsigned int x=0; x<V.cardinality(dim); x++)
            {
                if ((it=P.find(x))!=P.end())
                    res.push_back(it->second);
                else
                    res.push_back(0.0);

                //cout<<x<<":"<<*(res.end()-1)<<" ; ";
            }
            //cout<<")"<<endl;
            return res;
        }
        static int draw(const vector<double>& P, boost::mt19937& gen) {

            std::vector<double> cumulative;
            std::partial_sum(P.begin(), P.end(), std::back_inserter(cumulative));
//            cout<<"cumu = (";
//            for (unsigned int i=0; i<cumulative.size(); i++)
//                cout<<i<<":"<<cumulative[i]<<" ";
//            cout<<")"<<endl;
            boost::uniform_real<> dist(0, cumulative.back());
            boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(gen, dist);
            double draw=die();
    //        cout<<"unif draw = "<<draw<<endl;
            return (std::lower_bound(cumulative.begin(), cumulative.end(), draw) - cumulative.begin()) ;
        }

        static Distrib1D opposite(Distrib1D& P, Variable& V, unsigned int dim)
        {
            Distrib1D res;
            vector<double> P_vect(toVect(P, V, dim));
            for(unsigned int i=0; i<V.cardinality(dim); i++)
            {
                if((1-P_vect[i])>=FLT_EPSILON)
                    res[i]=(1-P_vect[i])/double(V.cardinality(dim)-1);
            }
            return res;

        }

    private:


};



class Simulator
{
  public:
  Simulator(Variable& state, Variable& action, double defaultDT_) : Xt(state), Ut(action), defaultDT(defaultDT_) {;}
    ContValue& simulate(double dt, unsigned int nb_step=1)
    {
        for(i=0; i<nb_step; i++)
            simulate(Xt, Ut, dt/double(nb_step));
        return Xt;
    }

    ContValue& simulate(ContValue& stateAction, unsigned int nb_step=1)
    {
        Xt.getValueVect().assign(stateAction.getValueVect().begin(), stateAction.getValueVect().begin()+Xt.getVar().nbDim());
        Ut.getValueVect().assign(stateAction.getValueVect().begin()+Xt.getVar().nbDim(), stateAction.getValueVect().end());
        return simulate(defaultDT, nb_step);
    }


    virtual void simulate(ContValue& state, ContValue& action, double dt)=0;

        /*
        //FIXME
    void test(string file_out,unsigned int t_max, double dt)
    {
        Xt[0]=0.01;
        Xt[1]=0.0;
        Ut[0]=0.0;

        vector<double> plotPos, plotSpeed, plotAction;
        for (unsigned int t=0; t<t_max; t++)
        {
            simulate(Xt,Ut,dt);
            plotPos.push_back(Xt[0]); plotSpeed.push_back(Xt[1]); plotAction.push_back(Ut[0]);

        }
        ofstream file(file_out.c_str());
        for (unsigned int t=0; t<t_max; t++)
        {
            file<<plotPos[t]<<" "<<plotSpeed[t]<<" "<<plotAction[t]<<" "<<endl;
        }
        file.close();
    }
        */

    ContValue& simulateGaussianDraw(double dt, vector<double>& sig, double nbSig, unsigned int nb_step=1)
    {
        Distributions::Distrib1D P;
        Value X_d(Xt.getVar());
        for(i=0; i<nb_step; i++)
            simulate(Xt, Ut, dt/double(nb_step));
      //  cout<<"Mean simulateGaussianDraw : ("<<Xt[0]<<" ; "<<Xt[1]<<")"<<endl;
        for (unsigned int i=0; i<Xt.getVar().nbDim(); i++)
        {
           // cout<<Xt[i]<<endl;
            Distributions::gaussian(P,Xt.getVar(),i,Xt[i],sig[i],nbSig);
            X_d[i]=Distributions::draw(Distributions::toVect(P,Xt.getVar(),i), gen);
            P.clear();
        }
        Xt.continuize(X_d);
        return Xt;

//        Distributions::Distrib1D P;
//        Value X_d(Xt.getVar());
//        simulate(Xt, Ut, dt, f);
//        Distributions::gaussian(P,Xt.getVar(),0,Xt[0],sigP,NB_SIG);
//        X_d[0]=Distributions::draw(Distributions::toVect(P,Xt.getVar(),0), gen);
//        P.clear();
//        Distributions::gaussian(P,Xt.getVar(),1,Xt[1],sigS,NB_SIG);
//        X_d[1]=Distributions::draw(Distributions::toVect(P,Xt.getVar(),1), gen);
//        Xt.continuize(X_d);
//        return Xt;

    }
    virtual ~Simulator();
    void setXU(ContValue& state, ContValue& action)
    {
        Xt=state;
        Ut=action;
    }
  protected:
    ContValue Xt,Ut;
    double defaultDT;
  private:
    boost::mt19937 gen;
    unsigned int i;
};


class ModelDescription
{
    public:
        ModelDescription(Variable& state_, Variable& action_, Simulator& sim_) : state(state_), action(action_), sim(sim_), cXt(state), cU(action), cXtp1(state), Xtp1(state) {;}
        void computeProbs();
        virtual double cost(const Value& state_, const Value& action_)=0;
        virtual double P_UkX(const Value& state_, const Value& action_)
        {
            return 1.0/double(action.cardinality());
        }

        virtual void getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U) =0;

        Variable& getStateVar() {return state;}
        Variable& getActionVar() {return action;}
        virtual ~ModelDescription();
    protected:
        Variable state,action;
        Simulator& sim;
        ContValue cXt;
        ContValue cU;
        ContValue cXtp1;
        Value Xtp1;
    private:
};






#endif // SIMULATOR_H
