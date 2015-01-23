#include "Learner.h"

void CndDistribution1DLearner::addPoint(ContValue& left, ContValue& right)
{
    leftVal_d.discretize(left);
    rightVal_d.discretize(right);

    if ( ! nb_data.count(rightVal_d))
        nb_data[rightVal_d]=init_weight;

    nb_data[rightVal_d]++;

    updateParams(left, right);




}

void CndHistogram1DLearner::updateParams(ContValue& left, ContValue& right)
{

    if (data[rightVal_d].find(leftVal_d)==data[rightVal_d].end()){

        data[rightVal_d][leftVal_d]=double(init_weight) * (sim.simulate(right))[leftDim];

    }
    data[rightVal_d][leftVal_d]++;

}

Distributions::Distrib1D CndHistogram1DLearner::getDistribution(const Value& rightValue)
{
    Distributions::Distrib1D res;
    std::map<Value, unsigned int>::iterator it;
    for(it=data[rightValue].begin(); it!=data[rightValue].end(); ++it)
    {
        res[(it->first)[leftDim]]=double(it->second)/double(nb_data[rightValue]);
    }
    return res;
}



void CndGaussian1DLearner::updateParams(ContValue& left, ContValue& right)
{
    unsigned int i;

    if (sum.find(rightVal_d)==sum.end()){
        sum[rightVal_d]= double(init_weight) * (sim.simulate(right))[leftDim];
        if (init_weight==0)
            squared_sum[rightVal_d]=0.0;
        else
            squared_sum[rightVal_d]= double(init_weight) * (init_sig*init_sig + (sum[rightVal_d]*sum[rightVal_d])/double(init_weight*init_weight));
    }
//    for(i=0; i<leftVar.nbDim(); i++)
//    {
        sum[rightVal_d]+=left[leftDim];
        squared_sum[rightVal_d]+=(left[leftDim]*left[leftDim]);

  //      cout<<"sum = "<<sum[rightVal_d][i]<< " ; sq_sum = "<<squared_sum[rightVal_d][i]<< " ; v ="<<left[i]<<" ; sq_v ="<<left[i]*left[i]<<"  ;  "<<endl;
//    }

}


Distributions::Distrib1D CndGaussian1DLearner::getDistribution(const Value& rightValue)
{
//    vector<Distributions::Distrib1D> gaussians(leftVar.nbDim());
    double mu, var;
    Distributions::Distrib1D res;

//    if (nb_data[rightValue]>=2)
//    {
//        for (unsigned int i=0; i<leftVar.nbDim(); i++)
//        {
            mu=sum[rightValue]/double(nb_data[rightValue]);
            var=squared_sum[rightValue]/double(nb_data[rightValue]) - mu*mu;
        //    cout<<"nb_data ="<<double(nb_data[rightValue])<<endl;
     //       cout<<"print : "<<sqrt(squared_sum[rightValue][i]/double(nb_data[rightValue]) - mu*mu)<<" ; "<<squared_sum[rightValue][i]<<" ; "<<double(nb_data[rightValue])<<endl;
            Distributions::gaussian(res, leftVar, leftDim, mu, sqrt(var), nbSig);
//        }
//
//        product(gaussians, gaussians.begin(), leftVal_d, res, 1.0);
//    }
           // cout<<"mu["<<leftDim<<"] = "<<mu<<endl;
            sum_sig+=sqrt(var);
            nb_dists++;
    return res;

}

void ModelLearner::addData(ContValue& left, ContValue& right)
{
    for(learners_it=learners.begin(); learners_it<learners.end(); learners_it++)
        (*learners_it)->addPoint(left, right);
}

void ModelLearner::update(Value& stateObj)
{
    Value actionVal(eng.action);
    Value stateVal(eng.state);
    vector<Distributions::Distrib1D> distribs1D;
    std::map<Value, unsigned int >::iterator right_it;

    for (right_it=learners[0]->nb_data.begin(); right_it!=learners[0]->nb_data.end(); ++right_it)   // to improve ...
    {
        for(learners_it=learners.begin(); learners_it<learners.end(); learners_it++)
            distribs1D.push_back((*learners_it)->getDistribution(right_it->first));
        Distributions::Distrib distrib;
        product(distribs1D, distribs1D.begin(), leftVal_d, distrib);
//        double sum=0.0;
//        for (Distributions::Distrib::iterator it=distrib.begin(); it!=distrib.end(); ++it)
//        {
//            cout<<"("<<it->first[0]<<","<<it->first[1]<<") : "<<it->second<<endl;
//            sum+=it->second;
//        }
//        cout<<"Sum = "<<sum<<" ; distrib.size() = "<<distrib.size()<<endl;
        distribs1D.clear();
//        if (distrib.size()!=0)
//        {
            stateVal.getValueVect().assign(right_it->first.getValueVect().begin(), right_it->first.getValueVect().begin()+eng.state.nbDim());
            actionVal.getValueVect().assign(right_it->first.getValueVect().begin()+eng.state.nbDim(), right_it->first.getValueVect().end());
            eng.removeDistrib(eng.getVertex(stateVal), actionVal);
            eng.insertDistribInGraph(distrib, eng.getVertex(stateVal), actionVal);
//        }


    }
    eng.optimize(stateObj);
    eng.toUpdate();
 //   eng.computePolicy(stateObj);

}

void ModelLearner::learnExhaustive(unsigned int nb_sample_per_right_val, unsigned int nb_step_sim, double dt)
{
    Value rval(rightVar);
    ContValue rval_c(rightVar);


    boost::mt19937 gen;
    gen.seed(time(NULL));
    boost::uniform_real<> unif(0, 1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > draw(gen, unif);

    vector<double> sig;
    sig.push_back(0.5);
    sig.push_back(0.5);
    rval.reset();
    unsigned int nb_right_processed=0;
    do
    {
        for (unsigned int i=0; i<nb_sample_per_right_val; i++)
        {
           // cout<<"rval_c=(";
            for(unsigned int dim=0; dim<rightVar.nbDim(); dim++)
            {
                rval_c[dim]=rightVar.getMins(dim)+(rval[dim]+draw())*(rightVar.getMaxs(dim)-rightVar.getMins(dim))/double(rightVar.cardinality(dim));
         //       cout<<rval_c[dim]<<", ";
            }
       //     cout<<endl;
//            sim.simulate(rval_c,0); // only to store state and action values in the simulator.
//            addData(sim.simulateGaussianDraw(dt,sig,4.0,nb_step_sim), rval_c);
            addData(sim.simulate(rval_c,nb_step_sim), rval_c);
        }
        cout<<"nb_right_processed = "<<(nb_right_processed++)<<endl;
    } while(rval.next());
}

void ModelLearnerOneGaussOneDirac::addData(ContValue& left, ContValue& right)
{
    x_c[0]=right[0];
    x_c[1]=left[1];

    //stateVal.getValueVect().assign(right.getValueVect().begin(), right.getValueVect().begin()+stateVal.getVar().nbDim());
    learners[0]->addPoint(left, x_c);
    learners[1]->addPoint(left, right);
}

void ModelLearnerOneGaussOneDirac::update(Value& stateObj)
{
    Value stateVal(eng.state);
    Value actionVal(eng.action);
    Value Xtp1(eng.state);
    Value PtVtp1(eng.state);
    ContValue Xtp1_c(eng.state);
    ContValue PtVtp1_c(eng.state);
    ContValue X_c(eng.state);
    std::map<Value, unsigned int >::iterator right_it;
    Distributions::Distrib1D::iterator pos_it, spd_it;

    for (right_it=learners[1]->nb_data.begin(); right_it!=learners[1]->nb_data.end(); ++right_it)
    {
        stateVal.getValueVect().assign(right_it->first.getValueVect().begin(), right_it->first.getValueVect().begin()+eng.state.nbDim());
        actionVal.getValueVect().assign(right_it->first.getValueVect().begin()+eng.state.nbDim(), right_it->first.getValueVect().end());
        PtVtp1[0]=stateVal[0];


        Distributions::Distrib1D spd=learners[1]->getDistribution(right_it->first);
        Distributions::Distrib distrib;
        for(spd_it=spd.begin(); spd_it!=spd.end(); ++spd_it)
        {
            PtVtp1[1]=(spd_it->first);

            Distributions::Distrib1D pos=learners[0]->getDistribution(PtVtp1);
            Xtp1[1]=PtVtp1[1];

/** version P(Ptp1 | Pt Vtp1) dirac :  **/
            PtVtp1_c.continuize(PtVtp1);
            Xtp1_c[0]=PtVtp1_c[0] + dt*PtVtp1_c[1];
            Xtp1_c[1]= PtVtp1_c[1];
            Xtp1.discretize(Xtp1_c);
            distrib[Xtp1]=spd_it->second;

/** version P(Ptp1 | Pt Vtp1) histogramme :  **/
//            for(pos_it=pos.begin(); pos_it!=pos.end(); ++pos_it)
//            {
//                Xtp1[0]=pos_it->first;
//
//                distrib[Xtp1]=pos_it->second * spd_it->second;
//            }

        }
        eng.removeDistrib(eng.getVertex(stateVal), actionVal);
        eng.insertDistribInGraph(distrib, eng.getVertex(stateVal), actionVal);
    }
//    cout<<"moy_sig = "<<dynamic_cast<CndGaussian1DLearner*>(learners[0])->sum_sig/double(dynamic_cast<CndGaussian1DLearner*>(learners[0])->nb_dists)<<endl;
//    dynamic_cast<CndGaussian1DLearner*>(learners[0])->sum_sig=0.0; dynamic_cast<CndGaussian1DLearner*>(learners[0])->nb_dists=0;
    eng.optimize(stateObj);
    eng.toUpdate();
}
