#ifndef LEARNER_H_INCLUDED
#define LEARNER_H_INCLUDED

#include "Optimizer.h"

class CndDistribution1DLearner
{
    public:
        CndDistribution1DLearner(Variable& leftVar_, unsigned int leftDim_, Variable& rightVar_, Simulator& sim_, unsigned int init_weight_) : leftVar(leftVar_), rightVar(rightVar_), leftVal_d(leftVar), rightVal_d(rightVar), sim(sim_), leftDim(leftDim_), init_weight(init_weight_) {;}
        void addPoint(ContValue& left, ContValue& right);
        virtual Distributions::Distrib1D getDistribution(const Value& rightValue)=0;

    protected:
    friend class ModelLearner;
    friend class ModelLearnerOneGaussOneDirac;
    Variable leftVar, rightVar;

    std::map<Value, unsigned int> nb_data;
    Value leftVal_d, rightVal_d;
    Simulator& sim;
    unsigned int leftDim;
    double init_weight;

    virtual void updateParams(ContValue& left, ContValue& right)=0;
};

class CndHistogram1DLearner : public CndDistribution1DLearner
{
    public:
        CndHistogram1DLearner(Variable leftVar_, unsigned int leftDim_, Variable rightVar_, Simulator& sim_, unsigned int init_weight_=0) : CndDistribution1DLearner(leftVar_, leftDim_, rightVar_, sim_,init_weight_) {;}
        Distributions::Distrib1D getDistribution(const Value& rightValue);
    private:
        std::map<Value, std::map<Value, unsigned int> > data;   /* associate right values to a map which associate to left value a counter */

        void updateParams(ContValue& left, ContValue& right);

};

class CndGaussian1DLearner : public CndDistribution1DLearner
{
    public:
        CndGaussian1DLearner(Variable leftVar_, unsigned int leftDim_, Variable rightVar_, double nbSig_, Simulator& sim_, double init_sig_, unsigned int init_weight_=0) : CndDistribution1DLearner(leftVar_, leftDim_, rightVar_, sim_, init_weight_), nbSig(nbSig_), init_sig(init_sig_), sum_sig(0.0), nb_dists(0) {;}
        Distributions::Distrib1D getDistribution(const Value& rightValue);
    private:
        friend class ModelLearnerOneGaussOneDirac;
        std::map<Value, double > sum;
        std::map<Value, double > squared_sum;
        double nbSig;
        double init_sig;
        void updateParams(ContValue& left, ContValue& right);
        double sum_sig;
        unsigned int nb_dists;



};


class ModelLearner
{
    public:
        ModelLearner(Optimizer& eng_, vector<CndDistribution1DLearner*> learners_, Variable& left, Variable right, Simulator& sim_) : eng(eng_), learners(learners_), leftVal_d(left), rightVar(right), sim(sim_) {;}
        virtual void addData(ContValue& left, ContValue& right);
        virtual void update(Value& stateObj);
        void learnExhaustive(unsigned int nb_sample_per_right_val, unsigned int nb_step_sim, double dt);
        //CndDistributionLearner& getLearner() {return learner;}
    protected:
        Optimizer& eng;
        vector<CndDistribution1DLearner*> learners;
        vector<CndDistribution1DLearner*>::iterator learners_it;
        Value leftVal_d;
        Variable rightVar;
        Simulator& sim;

        void product(vector<Distributions::Distrib1D>& dist1DVect, vector<Distributions::Distrib1D>::iterator it_dist1DVect, Value& v, Distributions::Distrib& dist, double p=1.0)
        {
//            cout<<"dist1DVect.size() = "<<dist1DVect.size()<<endl;
            for(Distributions::Distrib1D::iterator it_dist1D= it_dist1DVect->begin(); it_dist1D!=it_dist1DVect->end(); ++it_dist1D )
            {
                unsigned int indDistVect=it_dist1DVect-dist1DVect.begin();
//                cout<<"indDistVect = "<<indDistVect<<endl;
//                cout<<"it_dist1DVect->size() = "<<it_dist1DVect->size()<<endl;
                v[indDistVect]=it_dist1D->first;
                if (it_dist1DVect==dist1DVect.end()-1)
                    dist[v]=p*it_dist1D->second;
                else
                    product(dist1DVect, it_dist1DVect+1, v, dist, p*it_dist1D->second);
            }
        }

};


class ModelLearnerOneGaussOneDirac : public ModelLearner
{
    public:
    ModelLearnerOneGaussOneDirac(Optimizer& eng_, vector<CndDistribution1DLearner*> learners_, Variable& left, Variable right, Simulator& sim, double dt_) : ModelLearner(eng_,learners_,left, right, sim), x_c(eng_.getStateVar()), dt(dt_)
    {
        ;
    }
    void addData(ContValue& left, ContValue& right);
    void update(Value& stateObj);
    private:
        ContValue x_c;
        double dt;
};

//class TestLearner : public CndGaussian1DLearner
//{
//    public:
//        TestLearner(Variable leftVar_, Variable rightVar_, double nbSig_) : CndGaussian1DLearner(leftVar_, rightVar_, nbSig_, 0.0, 0), nbSig(nbSig_), leftVal_c(leftVar), rightVal_c(rightVar), file("learning.dat")
//        {
//            Distributions::gaussian(p0, leftVar, 0, 0.3, 0.1, nbSig);
//            Distributions::gaussian(p1, leftVar, 1, -0.4, 0.1, nbSig);
//        }
//
//        void learn(unsigned int nb_data)
//        {
//            rightVal_d[0]=0;
//            rightVal_c.continuize(rightVal_d);
//            for(unsigned int i=0; i<nb_data; i++)
//            {
//                leftVal_d[0]=Distributions::draw(Distributions::toVect(p0,leftVar,0), gen);
//                leftVal_d[1]=Distributions::draw(Distributions::toVect(p1,leftVar,1), gen);
//       //         cout<<"Discr : ("<<leftVal_d[0]<<","<<leftVal_d[1]<<") ; "<<endl;
//                leftVal_c.continuize(leftVal_d);
//       //         cout<<"Cont : ("<<leftVal_c[0]<<","<<leftVal_c[1]<<") ; "<<endl;
//                addPoint(leftVal_c, rightVal_c);
//            }
//            p_learn=getDistribution(rightVal_d);
//            leftVal_d.reset();
//            Distributions::Distrib::iterator it;
//            double sum=0.0;
//            do
//            {
//                if ((it=p_learn.find(leftVal_d))!=p_learn.end())
//                {
//                    file<<it->second<<" ";
//                    sum+=it->second;
//                }
//                else
//                    file<<0.0<<" ";
//            } while(leftVal_d.next());
////            for(Distributions::Distrib::iterator it=p_learn.begin(); it!=p_learn.end(); ++it)
////            {
////                file<<it->second<<" ";
////            }
//            file.close();
//            cout<<"SUM = "<<sum<<endl;
//
//        }
//    private:
//        double nbSig;
//        Distributions::Distrib1D p0;
//        Distributions::Distrib1D p1;
//        Distributions::Distrib p_learn;
//        ContValue leftVal_c, rightVal_c;
//        boost::mt19937 gen;
//        ofstream file;
//};

#endif // LEARNER_H_INCLUDED
