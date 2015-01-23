#include "Leg2DHDescription.h"


void Leg2DHDescription::getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U)
{
    Distributions::Distrib1D p_theta, p_dtheta,p_h, p_dh;

    cXt.continuize(Xt);
    cU.continuize(U);
    sim.setXU(cXt,cU);
    cXtp1=sim.simulate(DT);

//    vector<double> mus_sigs(learnParams(Xt,U,1));

    Distributions::gaussian(p_theta, state, 0, cXtp1[0], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_dtheta, state, 1, cXtp1[1], SIGMA_V, NB_SIG);
    Distributions::gaussian(p_h, state, 2, cXtp1[2], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_dh, state, 3, cXtp1[3], SIGMA_V, NB_SIG);
    
    Distributions::Distrib1D::iterator theta_it, dtheta_it,h_it, dh_it;
    Value Xtp1(state);

    
    
    for(theta_it=p_theta.begin(); theta_it!=p_theta.end(); ++theta_it)
    {
        for(dtheta_it=p_dtheta.begin(); dtheta_it!=p_dtheta.end(); ++dtheta_it)
        {   
            for(h_it=p_h.begin(); h_it!=p_h.end(); ++h_it)
            {

                for(dh_it=p_dh.begin(); dh_it!=p_dh.end(); ++dh_it)
                {                
            
                    Xtp1[0]=theta_it->first;
                    Xtp1[1]=dtheta_it->first;
                    Xtp1[2]=h_it->first;
                    Xtp1[3]=dh_it->first;
                    
                    P[Xtp1]=theta_it->second * dtheta_it->second * h_it->second * dh_it->second;
                }
            }
            
        }
    }
        /*
    if(xs<=FLT_EPSILON || ys<=FLT_EPSILON ||thetas<=FLT_EPSILON)
        cerr<<"DISTRIB ERROR "<<Xt[0]<<" "<<Xt[1]<<" "<<Xt[2]<<" "<<Xt[3]<<" "<<U[0]<<endl;
        */
}


Leg2DHDescription::~Leg2DHDescription()
{
        //dtor
}

void Leg2DHOneGaussDescription::getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U)
{
    Distributions::Distrib1D p_theta, p_dtheta,p_h, p_dh;
    cXt.continuize(Xt);
    cU.continuize(U);
    sim.setXU(cXt,cU);
    cXtp1=sim.simulate(DT,NB_STEP_SIM);

    Distributions::gaussian(p_theta, state, 0, cXtp1[0], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_dtheta, state, 1, cXtp1[1], SIGMA_V, NB_SIG);
    Distributions::gaussian(p_h, state, 2, cXtp1[2], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_dh, state, 3, cXtp1[3], SIGMA_V, NB_SIG);
    
    Distributions::Distrib1D::iterator theta_it, dtheta_it,h_it, dh_it;
    Value Xtp1(state);

    
    
    for(theta_it=p_theta.begin(); theta_it!=p_theta.end(); ++theta_it)
    {
        for(dtheta_it=p_dtheta.begin(); dtheta_it!=p_dtheta.end(); ++dtheta_it)
        {   
            for(h_it=p_h.begin(); h_it!=p_h.end(); ++h_it)
            {

                for(dh_it=p_dh.begin(); dh_it!=p_dh.end(); ++dh_it)
                {                
            
                    Xtp1[0]=theta_it->first;
                    Xtp1[1]=dtheta_it->first;
                    Xtp1[2]=h_it->first;
                    Xtp1[3]=dh_it->first;
                    
                    P[Xtp1]=theta_it->second * dtheta_it->second * h_it->second * dh_it->second;
                }
            }
            
        }
    }
    
}

void Leg2DHSimulator::simulate(ContValue& X, ContValue& U, double dt)
{


        /* Inverse kinematics
         * H_dot=J.Q_dot
         *
         * H=|theta|
         *   |  h  |
         *
         * SIMPLE version controle en vitesse (couple infini). dtheta ne dépend pas de dtheta_t-1
         * 
         * Q_dot=|alpha|
         *       |beta |
         *
         * J=|dtheta/dalpha  dtheta/dbeta|
         *   |dh/dalpha      dh/dbeta|
         * 
         */
    
    double h=X[2];
    double theta=X[0];
    
    
    double beta=0.0;
        //double beta=acos((L1*L1+L2*L2-h*h)/(2.0*L1*L2));    
    
    if(h>=L1+L2)
        beta=M_PI;
    else
        beta=acos((L1*L1+L2*L2-h*h)/(2.0*L1*L2));    


    double alpha=M_PI/2.0-beta/2.0+atan((L2-L1)/(L2+L1)*cotan(beta))+theta;

        //1 dof
        //double dalpha=U[0];
    double dbeta=U[0];
    
    
        
        //double J11=1.0;
    
    double J12=-(-L1 + L2)*(-0.5*cotan(0.5*beta)*cotan(0.5*beta) - 0.5)/((L1 + L2)*((-L1 + L2)*(-L1 + L2)*cotan(0.5*beta)*cotan(0.5*beta)/(L1 + L2)*(L1 + L2) + 1)) + 0.5;
    

        //double J21=0;
    double J22=L1*L2*sin(beta)/sqrt(L1*L1 - 2.0*L1*L2*cos(beta) + L2*L2);


    double prev_dbeta=X[3]/J22;
    double prev_dalpha=X[1]-prev_dbeta*J12;
    
    double dht=dbeta*J22;
//    double dbeta=dht/J22;
    
        //constant dalpha
    double dthetat=prev_dalpha+dbeta*J12;

    double thetat=theta+dthetat*dt;
    double ht=h+dht*dt;

        /*
        //Si pied au sol, pas de mouvement
    if(ht>=L1+L2-MINHFOOT){
        dthetat=0.0;
        thetat=theta;
        
    }
        */
    if(ht>=HMAX){    
        ht=HMAX-EPSILON;
        if(dht>0.0)
            dht=0.0;
        
    }

    if(ht<=HMIN){    
        ht=HMIN+EPSILON;
        if(dht<0.0)
            dht=0.0;
        
    }

    if(thetat<THETAMIN){    
        thetat=THETAMIN;
        if(dthetat<0.0)
            dthetat=0.0;
        
    }
    if(thetat>THETAMAX){    
        thetat=THETAMAX;
        if(dthetat>0.0)
            dthetat=0.0;
        
    }

    
    
        //double beta=arccos((L2_2+L1_1-ht*ht)/(2.0*L1*L2)  );        

        /*
    if(h==1.99)
        printf("DEBUG: %f %f %f %f %f %f %f %f\n",alpha,beta,prev_dalpha,dbeta,thetat,dthetat,ht,dht);
        */
    
    X[0]=thetat;
    X[1]=dthetat;
    X[2]=ht;
    X[3]=dht;
    
    for (unsigned int i=0; i<X.getVar().nbDim(); i++)
        X.circularize(i);
}


void Leg2DHManager::buildVariables()
{
    vector<unsigned int> stateCards(4);
    vector<double> stateMins(4);
    vector<double> stateMaxs(4);
    vector<bool> stateCircular(4);
        //stateCards[0]=XCARD; stateMins[0]=XMIN; stateMaxs[0]=XMAX; stateCircular[0]=X_CIRC;
    stateCards[0]=THETACARD; stateMins[0]=THETAMIN; stateMaxs[0]=THETAMAX; stateCircular[0]=false;
    stateCards[1]=DTHETACARD; stateMins[1]=DTHETAMIN; stateMaxs[1]=DTHETAMAX; stateCircular[1]=false;
    stateCards[2]=HCARD; stateMins[2]=HMIN; stateMaxs[2]=HMAX; stateCircular[2]=false;
    stateCards[3]=DHCARD; stateMins[3]=DHMIN; stateMaxs[3]=DHMAX; stateCircular[3]=false;

    
    state = new Variable(stateMins, stateMaxs, stateCards, stateCircular);

        /*
    vector<unsigned int> actionCards(2);
    vector<double> actionMins(2);
    vector<double> actionMaxs(2);
    vector<bool> actionCircular(2);
    actionCards[0]=DALPHACARD; actionMins[0]=DALPHAMIN; actionMaxs[0]=DALPHAMAX; actionCircular[0]=false;
    actionCards[1]=DBETACARD; actionMins[1]=DBETAMIN; actionMaxs[1]=DBETAMAX; actionCircular[0]=false;
        */

    vector<unsigned int> actionCards(1);
    vector<double> actionMins(1);
    vector<double> actionMaxs(1);
    vector<bool> actionCircular(1);
    actionCards[0]=DALPHACARD; actionMins[0]=DALPHAMIN; actionMaxs[0]=DALPHAMAX; actionCircular[0]=false;
    action= new Variable(actionMins, actionMaxs, actionCards, actionCircular);


}

void Leg2DHManager::buildSimulator()
{
        //sim = new Leg2DHSimulator(*state, *action, DT, XMIN, XMAX);
    sim = new Leg2DHSimulator(*state, *action, DT);
}

void Leg2DHManager::buildModel()
{
    model = new Leg2DHDescription(*state, *action, *sim);
}

void Leg2DHManager::buildOptimizer(bool drawnPolicy, double initBeta)
{
    optimizer = new OptimizerQuasiDist(*model, drawnPolicy, initBeta);
}

void Leg2DHManager::buildDistribLearner()
{
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 0, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 1, (*state)^(*action), NB_SIG, *sim, SIGMA_V, LEARN_INIT_WEIGHT));
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 2, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 3, (*state)^(*action), NB_SIG, *sim, SIGMA_V, LEARN_INIT_WEIGHT));
}


void Leg2DHManager::buildModelLearner()
{
    learner = new ModelLearner(*optimizer, distrib1DLearners, *state, (*state)^(*action), *sim);
}

void Leg2DHManager::buildController()
{
    controller = new Controller(*optimizer, *learner);
}

void Leg2DHManager::buildObjective()
{
    (*obj_c)[0]=OBJ_THETA;
    (*obj_c)[1]=OBJ_DTHETA;
    (*obj_c)[2]=OBJ_H;
    (*obj_c)[3]=OBJ_DH;
    
}

void Leg2DHManager::run()
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
    vector<double> plotTheta, plotDTheta, plotH, plotDH, plotDAlpha,plotDBeta;
    vector<unsigned int> plotThetaD, plotDThetaD, plotHD, plotDHD, plotDAlphaD,plotDBetaD;

    vector<double> sigmas;  /* standard deviation during the control simulation below */
    sigmas.push_back(0.0);

    sigmas.push_back(0.0);

    sigmas.push_back(0.0);
    sigmas.push_back(0.0);
    
    
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


        //int simmap[51][51]={0};
    
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


        /*
    ofstream fpolicy("full_policy.dat");
    
    for(int i=0;i<THETACARD;i++)
    {
        for(int j=0;j<DTHETACARD;j++)
        {
            for(int k=0;k<HCARD;k++)
            {
                for(int l=0;l<DHCARD;l++)
                {
                
                Xtmp[0]=i;
                Xtmp[1]=j;
                Xtmp[2]=k;
                Xtmp[3]=l;
                
                
                fpolicy<<i<<" "<<j<<" "<<k<<" "<<l<<" ";
                
                for(int m=0;m<DALPHACARD;m++)
                {
                    for(int n=0;m<DBETACARD;n++)
                    {
                    
                    fpolicy<<optimizer->getPolicy(optimizer->getVertex(Xtmp))[m][n]<<" ";

                }
                fpolicy<<endl;

            }
            
        }
        
    }
    
    
    fpolicy.close();
    
        */




    
    for (unsigned int nb_sim=0; nb_sim<1; nb_sim++)
    {

        cout<<"sim: "<<nb_sim<<endl;
        
    
        Xt_c[0]=0.0;
        Xt_c[1]=0.0;
        Xt_c[2]=L1+L2;
        Xt_c[3]=0.0;
   
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

                    //simmap[Xt[0]][Xt[1]]+=1;
                
                

                U=controller->control(Xt, *obj);
                
                    

                plotThetaD.push_back(Xt[0]);
                plotDThetaD.push_back(Xt[1]);
                plotHD.push_back(Xt[2]);
                plotDHD.push_back(Xt[3]);
                plotDBetaD.push_back(U[0]);
                    //plotDBetaD.push_back(U[1]);
                
                U_c.continuize(U);
       //         cout<<t<<" ; ("<<Xt[0]<<", "<<Xt[1]<<", "<<U[0]<<")  ;  ("<<Xt_c[0]<<", "<<Xt_c[1]<<", "<<U_c[0]<<")"<<endl;


                plotTheta.push_back(Xt_c[0]);
                plotDTheta.push_back(Xt_c[1]);
                plotH.push_back(Xt_c[2]);
                plotDH.push_back(Xt_c[3]);
                plotDBeta.push_back(U_c[0]);
                    //plotDBeta.push_back(U_c[1]);
                
                    

               // dynamic_cast<Leg2DHSimulator*>(sim)->setFriction(0.2);
          //      dynamic_cast<Leg2DHSimulator*>(sim)->setSlope(0.18);



                sim->setXU(Xt_c, U_c);


           /***     Deterministic sim ****/
                Xtp1_c=sim->simulate(DT,NB_STEP_SIM);
//                cout<<"deterministic : "<<Xtp1_c[0]<<", "<<Xtp1_c[1]<<endl;

           /***     Gaussian noised sim ****/
                    //Xtp1_c=sim->simulateGaussianDraw(DT, sigmas,5.0, NB_STEP_SIM);
                

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
            file<<plotTheta[t]<<" "<<plotDTheta[t]<<" "<<plotH[t]<<" "<<plotDH[t]<<" "<<plotDBeta[t]/*<<" "<<plotDBeta[t]*/<<endl;
        }
        file.close();

        ofstream file_discr((folder+"/control_discr.dat").c_str());
        for (unsigned int t=0; t<max_t_plot; t++)
        {
            file_discr<<plotThetaD[t]<<" "<<plotDThetaD[t]<<" "<<plotHD[t]<<" "<<plotDHD[t]<<" "<<plotDBetaD[t]/*<<" "<<plotDBetaD[t]*/<<endl;
        }

        file_discr.close();

                    
        
           /***     landscape and policy saving for plot ****/
        optimizer->toPlot(*obj);




}

           /***     Parameter learning from previously stored data (if needed) ****///        learner->update(*obj);
//        learner->update(*obj);
//        optimizer->computePolicy(*obj);
//        optimizer->toPlot(*obj);



Leg2DHManager::~Leg2DHManager()
{
    //dtor
}
