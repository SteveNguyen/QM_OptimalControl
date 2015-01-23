#include "NonHoloDescription.h"


void NonHoloDescription::getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U)
{
    Distributions::Distrib1D p_xpos, p_ypos,p_theta;

    cXt.continuize(Xt);
    cU.continuize(U);
    sim.setXU(cXt,cU);
    cXtp1=sim.simulate(DT);

//    vector<double> mus_sigs(learnParams(Xt,U,1));

    Distributions::gaussian(p_xpos, state, 0, cXtp1[0], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_ypos, state, 1, cXtp1[1], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_theta, state, 2, cXtp1[2], SIGMA_P, NB_SIG);
    
    Distributions::Distrib1D::iterator xpos_it, ypos_it,theta_it;
    Value Xtp1(state);

    float xs=0.0,ys=0.0,thetas=0.0;
    
    
    for(xpos_it=p_xpos.begin(); xpos_it!=p_xpos.end(); ++xpos_it)
    {
        xs+=xpos_it->second;
        
        for(ypos_it=p_ypos.begin(); ypos_it!=p_ypos.end(); ++ypos_it)
        {
            ys+=ypos_it->second;
            
            for(theta_it=p_theta.begin(); theta_it!=p_theta.end(); ++theta_it)
            {
                thetas+=theta_it->second;
                
            
                Xtp1[0]=xpos_it->first;
                Xtp1[1]=ypos_it->first;
                Xtp1[2]=theta_it->first;
                P[Xtp1]=xpos_it->second * ypos_it->second * theta_it->second;
                    //std::cerr<<Xt[0]<<" "<<Xt[1]<<" "<<U[0]<<" "<<Xtp1[0]<<" "<<Xtp1[1]<<" "<<P[Xtp1]<<std::endl;

                //std::cerr<<Xt[0]<<" "<<Xt[1]<<" "<<U[0]<<" "<<Xtp1[0]<<" "<<Xtp1[1]<<" "<<pos_it->second<<" "<<spd_it->second<<std::endl;
            }
        }
    }
    if(xs<=FLT_EPSILON || ys<=FLT_EPSILON ||thetas<=FLT_EPSILON)
        cerr<<"DISTRIB ERROR "<<Xt[0]<<" "<<Xt[1]<<" "<<Xt[2]<<" "<<U[0]<<endl;
    
}


//vector<double> NonHoloDescription::learnParams(Value& Xt, Value& U, unsigned int nb_sample_per_dim)
//{
//    cXt.continuize(Xt);
//    cU.continuize(U);
//    double sum_x0=0.0; double square_sum_x0=0.0;
//    double sum_cos_x0=0.0; double sum_sin_x0=0.0;
//    double sum_x1=0.0; double square_sum_x1=0.0;
//    for (unsigned int i_u=0; i_u<nb_sample_per_dim; i_u++)
//    {
//        cU[0]+=(double(i_u)-(double(nb_sample_per_dim)/2.0-0.5))*((action.getMaxs(0)-action.getMins(0))/double(nb_sample_per_dim*action.cardinality(0)));
//        for (unsigned int i_x0=0; i_x0<nb_sample_per_dim; i_x0++)
//        {
//            cXt[0]+=(double(i_x0)-(double(nb_sample_per_dim)/2.0-0.5))*((state.getMaxs(0)-state.getMins(0))/double(nb_sample_per_dim*state.cardinality(0)));
//            for (unsigned int i_x1=0; i_x1<nb_sample_per_dim; i_x1++)
//            {
//                cXt[1]+=(double(i_x1)-(double(nb_sample_per_dim)/2.0-0.5))*((state.getMaxs(1)-state.getMins(1))/double(nb_sample_per_dim*state.cardinality(1)));
//                sim.setXU(cXt,cU);
//                cXtp1=sim.simulate(DT);
//                sum_sin_x0+=sin(cXtp1[0]);
//                sum_cos_x0+=cos(cXtp1[0]);
//               // square_sum_x0+=cXtp1[0]*cXtp1[0];
//                sum_x1+=cXtp1[1];
//                square_sum_x1+=cXtp1[1]*cXtp1[1];
//
//            }
//
//
//        }
//    }
//
//    double mean_x0=atan2(sum_sin_x0, sum_cos_x0);
////    if (mean_x0<0.0)  // now useless due to the new angle range ( = [-pi, pi], ie the same than atan2)
////        mean_x0+=2.0*M_PI;
//    double mean_resultant_length=sqrt(sum_cos_x0*sum_cos_x0 + sum_sin_x0*sum_sin_x0)/pow(nb_sample_per_dim, state.nbDim()+action.nbDim());
//    double conc;
//    if (mean_resultant_length<0.53)
//        conc=2.0*mean_resultant_length+pow(mean_resultant_length,3)+5.0*pow(mean_resultant_length,5)/6.0;
//    else if (mean_resultant_length<0.85)
//        conc=-0.4 + 1.39*mean_resultant_length + 0.43/(1.0-mean_resultant_length);
//    else
//        conc=1.0/(pow(mean_resultant_length,3) - 4.0*pow(mean_resultant_length,2) + 3.0*mean_resultant_length);
//    conc=10.0;
// //   double std_x0=sqrt(square_sum_x0/double(nb_sample*nb_sample*nb_sample)-mean_x0*mean_x0);
//    double mean_x1=sum_x1/pow(nb_sample_per_dim, state.nbDim()+action.nbDim());
//    double std_x1=sqrt(square_sum_x1/pow(nb_sample_per_dim, state.nbDim()+action.nbDim())- pow(mean_x1,2));
//  //  cout<<"("<<U[0]<<","<<Xt[0]<<","<<Xt[1]<<"). Mean = ("<<mean_resultant_length<<","<<mean_x1<<") ; Std = ("<<conc<<","<<std_x1<<")."<<endl;
//
//  vector<double> res;
//  res.push_back(mean_x0); res.push_back(mean_x1); res.push_back(conc); res.push_back(std_x1);
//  return res;
//}

NonHoloDescription::~NonHoloDescription()
{
        //dtor
}

void NonHoloOneGaussDescription::getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U)
{
    cXt.continuize(Xt);
    cU.continuize(U);
    sim.setXU(cXt,cU);
    cXtp1=sim.simulate(DT,NB_STEP_SIM);

    Distributions::Distrib1D p_xpos, p_ypos, p_theta;
        //Distributions::dirac(p_pos, state, 0, cXt[0]+DT*cXtp1[1]);

    Distributions::gaussian(p_xpos, state, 0, cXtp1[0], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_ypos, state, 1, cXtp1[1], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_theta, state, 2, cXtp1[2], SIGMA_P, NB_SIG);
    

    Distributions::Distrib1D::iterator xpos_it, ypos_it,theta_it;
    Value Xtp1(state);

    for(xpos_it=p_xpos.begin(); xpos_it!=p_xpos.end(); ++xpos_it)
    {
        for(ypos_it=p_ypos.begin(); ypos_it!=p_ypos.end(); ++ypos_it)
        {
            for(theta_it=p_theta.begin(); theta_it!=p_theta.end(); ++theta_it)
            {
                
            
                Xtp1[0]=xpos_it->first;
                Xtp1[1]=ypos_it->first;
                Xtp1[2]=theta_it->first;
                P[Xtp1]=xpos_it->second * ypos_it->second * theta_it->second;
                
            }
        }
    }
    
}

void NonHoloSimulator::simulate(ContValue& X, ContValue& U, double dt)
{
        /*
    xt=X[0]+cos(X[2])*U[0]*DT
    yt=X[1]+sin(X[2])*U[0]*DT
    phi=X[2]+U[1]*DT

    if abs(phi)>pi:
        if phi<0.0:
            phi=phi+2.0*pi
        else:
            phi=phi-2.0*pi


    return (xt,yt,phi)
        */


        /*
    double xt=X[0]+cos(X[2])*dt;
    double yt=X[1]+sin(X[2])*dt;
    double thetat=X[2]+U[0]*dt;
        */        
        

        
    double xt=X[0];
    double yt=X[1];
    double thetat=X[2];

    for(int i=0;i<100;i++)
    {
        xt+=cos(thetat)*dt/100.0;
        yt+=sin(thetat)*dt/100.0;
        thetat+=U[0]*dt/100.0;
    }
        

    
        /*
    if(thetat>M_PI)
        thetat=thetat-2.0*M_PI;
    else if(thetat<=-M_PI)
        thetat=2.0*M_PI+thetat;
        */

        /*
    if(thetat>2.0*M_PI)
        thetat=thetat-2.0*M_PI;
    else if(thetat<=0.0)
        thetat=2.0*M_PI+thetat;
        */
    
    X[0]=xt;
    X[1]=yt;
    X[2]=thetat;
    
    for (unsigned int i=0; i<X.getVar().nbDim(); i++)
        X.circularize(i);
}


void NonHoloManager::buildVariables()
{
    vector<unsigned int> stateCards(3);
    vector<double> stateMins(3);
    vector<double> stateMaxs(3);
    vector<bool> stateCircular(3);
        //stateCards[0]=XCARD; stateMins[0]=XMIN; stateMaxs[0]=XMAX; stateCircular[0]=X_CIRC;
    stateCards[0]=XCARD; stateMins[0]=XMIN; stateMaxs[0]=XMAX; stateCircular[0]=false;
    stateCards[1]=YCARD; stateMins[1]=YMIN; stateMaxs[1]=YMAX; stateCircular[1]=false;
    stateCards[2]=THETACARD; stateMins[2]=THETAMIN; stateMaxs[2]=THETAMAX; stateCircular[2]=true;

    
    state = new Variable(stateMins, stateMaxs, stateCards, stateCircular);

    vector<unsigned int> actionCards(1);
    vector<double> actionMins(1);
    vector<double> actionMaxs(1);
    vector<bool> actionCircular(1);
    actionCards[0]=UCARD; actionMins[0]=UMIN; actionMaxs[0]=UMAX; actionCircular[0]=false;

    action= new Variable(actionMins, actionMaxs, actionCards, actionCircular);


}

void NonHoloManager::buildSimulator()
{
        //sim = new NonHoloSimulator(*state, *action, DT, XMIN, XMAX);
    sim = new NonHoloSimulator(*state, *action, DT);
}

void NonHoloManager::buildModel()
{
    model = new NonHoloDescription(*state, *action, *sim);
}

void NonHoloManager::buildOptimizer(bool drawnPolicy, double initBeta)
{
    optimizer = new OptimizerQuasiDist(*model, drawnPolicy, initBeta);
}

void NonHoloManager::buildDistribLearner()
{
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 0, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 1, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 2, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
}

void NonHoloManager::buildModelLearner()
{
    learner = new ModelLearner(*optimizer, distrib1DLearners, *state, (*state)^(*action), *sim);
}

void NonHoloManager::buildController()
{
    controller = new Controller(*optimizer, *learner);
}

void NonHoloManager::buildObjective()
{
    (*obj_c)[0]=OBJ_PX;
    (*obj_c)[1]=OBJ_PY;
    (*obj_c)[2]=OBJ_THETA;
    
}

void NonHoloManager::run()
{
//    cout<<"Learning ..."<<endl;
//    learner->learnExhaustive(100,100, DT);
//    learner->update(*obj);
//    cout<<"Done."<<endl;
//    optimizer->computePolicy(*obj);
//    optimizer->toPlot(*obj);
//    save("ex_learn.ser");
  //  exit(0);

    unsigned int max_t=500;
    unsigned int max_t_plot=500;
    Value Xt(*state);
    Value Xtp1(*state);
    Value U(*action);
    ContValue U_c(*action);
    ContValue Xt_c(*state);
    ContValue Xtp1_c(*state);
    Variable stateAction((*state)^(*action));

    ContValue XtU_c(stateAction);
    vector<double> plotX, plotY, plotTheta, plotAction;
    vector<unsigned int> plotXDiscr, plotYDiscr, plotThetaDiscr,plotActionDiscr;

    vector<double> sigmas;  /* standard deviation during the control simulation below */
    sigmas.push_back(0.0) /* x position*/;
        //sigmas.push_back(0.1); /* speed */
    sigmas.push_back(0.0); /* y position */

    sigmas.push_back(0.0); /* theta */
    
           /***     Folder where to store logs ****/
    string folder(".");
    optimizer->setDir(folder);

//    optimizer->toPlot(*obj);

//    optimizer->nbEdgesModif=0;
//    unsigned int nb_runs=1;

        //Xt_c[0]=-M_PI;
        //Xt_c[1]=0.0;

        //Xt_c[0]=M_PI/2.0;
        //Xt_c[1]=1.875;


    int simmap[51][51]={0};
    
        /*

    Xt_c[0]=0.0;
    Xt_c[1]=0.0;
    Xt_c[2]=0.0;
        //Xt_c[2]=-M_PI;
    
    Xt.discretize(Xt_c);
        */


//    optimizer->incrBeta(261);
//    cout<<"BETA = "<<optimizer->getBeta()<<endl;
//    optimizer->computePolicy(*obj);



        //Dump the full policy

    Value Xtmp(*state);
    ofstream fpolicy("full_policy.dat");
    
    for(int i=0;i<XCARD;i++)
    {
        for(int j=0;j<YCARD;j++)
        {
            for(int k=0;k<THETACARD;k++)
            {
                Xtmp[0]=i;
                Xtmp[1]=j;
                Xtmp[2]=k;
                
                fpolicy<<i<<" "<<j<<" "<<k<<" ";
                
                for(int l=0;l<UCARD;l++)
                {
                    
                    fpolicy<<optimizer->getPolicy(optimizer->getVertex(Xtmp))[l]<<" ";

                }
                fpolicy<<endl;

            }
            
        }
        
    }
    
    
    fpolicy.close();
    





    
    for (unsigned int nb_sim=0; nb_sim<1000; nb_sim++)
    {

        cout<<"sim: "<<nb_sim<<endl;
        
    
        Xt_c[0]=0.0;
        Xt_c[1]=0.0;
        Xt_c[2]=0.0;
            //Xt_c[2]=-M_PI;
        
        Xt.discretize(Xt_c);
        
        
        
//        if (nb_sim>0)
//        {
//            optimizer->nbEdgesModif=0;
//            folder="with_learning";
//            optimizer->setDir(folder);
//            //optimizer->toPlot(*obj);
//            nb_runs=1;
//
//        }
//
//    Xt_c[0]=-M_PI; //0.0;
//    Xt_c[1]=0.0;//(3.0/4.0)*VMAX;
//    Xt.discretize(Xt_c);
//
//        for (unsigned int nb_learn_runs=0; nb_learn_runs<nb_runs; nb_learn_runs++)
//        {
//
//
//            plotPos.clear(); plotSpeed.clear(); plotAction.clear();
//            plotPosDiscr.clear(); plotSpeedDiscr.clear(); plotActionDiscr.clear();
            for (unsigned int t=0; t<max_t; t++)
            {
            //    cout<<"t="<<t<<" ; U="<<U[0]<<" ; Uc="<<U_c[0]<<endl;


                    /*

	      /*

                for(int i=0;i<11;i++)
                cout<<optimizer->getPolicy(optimizer->getVertex(Xt))[i]<<" ";
                cout<<endl;

                    */

                simmap[Xt[0]][Xt[1]]+=1;
                
                

                U=controller->control(Xt, *obj);
                
                    

                plotXDiscr.push_back(Xt[0]);
                plotYDiscr.push_back(Xt[1]);
                plotThetaDiscr.push_back(Xt[2]);
                plotActionDiscr.push_back(U[0]);
                U_c.continuize(U);
       //         cout<<t<<" ; ("<<Xt[0]<<", "<<Xt[1]<<", "<<U[0]<<")  ;  ("<<Xt_c[0]<<", "<<Xt_c[1]<<", "<<U_c[0]<<")"<<endl;
                plotX.push_back(Xt_c[0]);
                plotY.push_back(Xt_c[1]);
                plotTheta.push_back(Xt_c[2]);
                plotAction.push_back(U_c[0]);

                    

               // dynamic_cast<NonHoloSimulator*>(sim)->setFriction(0.2);
          //      dynamic_cast<NonHoloSimulator*>(sim)->setSlope(0.18);



                sim->setXU(Xt_c, U_c);


           /***     Deterministic sim ****/
                    //Xtp1_c=sim->simulate(DT,NB_STEP_SIM);
//                cout<<"deterministic : "<<Xtp1_c[0]<<", "<<Xtp1_c[1]<<endl;

           /***     Gaussian noised sim ****/
                Xtp1_c=sim->simulateGaussianDraw(DT, sigmas,5.0, NB_STEP_SIM);


                Xtp1.discretize(Xtp1_c);


           /***     Adding data to learn (if needed) ****/
//                XtU_c=Xt_c^U_c;
//                learner->addData(Xtp1_c, XtU_c);


                Xt_c=Xtp1_c;
                Xt=Xtp1;




            }
            
//
        }


           /***     Sim logs ****/


                    
        ofstream file((folder+"/control.dat").c_str());
        for (unsigned int t=0; t<max_t_plot; t++)
        {
            file<<plotX[t]<<" "<<plotY[t]<<" "<<plotTheta[t]<<" "<<plotAction[t]<<" "<<endl;
        }
        file.close();

        ofstream file_discr((folder+"/control_discr.dat").c_str());
        for (unsigned int t=0; t<max_t_plot; t++)
        {
            file_discr<<plotXDiscr[t]<<" "<<plotYDiscr[t]<<" "<<plotThetaDiscr[t]<<" "<<plotActionDiscr[t]<<" "<<endl;
        }

        file_discr.close();

                    
        
           /***     landscape and policy saving for plot ****/
        optimizer->toPlot(*obj);




}

           /***     Parameter learning from previously stored data (if needed) ****///        learner->update(*obj);
//        learner->update(*obj);
//        optimizer->computePolicy(*obj);
//        optimizer->toPlot(*obj);



NonHoloManager::~NonHoloManager()
{
    //dtor
}
