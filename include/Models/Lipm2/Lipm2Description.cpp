#include "Lipm2Description.h"


void Lipm2Description::getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U)
{
    Distributions::Distrib1D p_pos, p_spd, f_pos, f_spd;

    cXt.continuize(Xt);
    cU.continuize(U);
    sim.setXU(cXt,cU);
    cXtp1=sim.simulate(DT);

//    vector<double> mus_sigs(learnParams(Xt,U,1));

    Distributions::gaussian(p_pos, state, 0, cXtp1[0], SIGMA_P, NB_SIG);
    Distributions::gaussian(p_spd, state, 1, cXtp1[1], SIGMA_S, NB_SIG);
    Distributions::gaussian(f_pos, state, 2, cXtp1[2], SIGMA_P, NB_SIG);
    Distributions::gaussian(f_spd, state, 3, cXtp1[3], SIGMA_S, NB_SIG);
    
    Distributions::Distrib1D::iterator pos_it, spd_it, fpos_it, fspd_it;
    Value Xtp1(state);
    for(pos_it=p_pos.begin(); pos_it!=p_pos.end(); ++pos_it)
    {
        for(spd_it=p_spd.begin(); spd_it!=p_spd.end(); ++spd_it)
        {

            for(fpos_it=f_pos.begin(); fpos_it!=f_pos.end(); ++fpos_it)
            {
                for(fspd_it=f_spd.begin(); fspd_it!=f_spd.end(); ++fspd_it)
                {

                    
                    Xtp1[0]=pos_it->first;
                    Xtp1[1]=spd_it->first;
                    Xtp1[2]=fpos_it->first;
                    Xtp1[3]=fspd_it->first;
                    
                    P[Xtp1]=pos_it->second * spd_it->second*fpos_it->second * fspd_it->second;
                //std::cerr<<Xt[0]<<" "<<Xt[1]<<" "<<U[0]<<" "<<Xtp1[0]<<" "<<Xtp1[1]<<" "<<P[Xtp1]<<std::endl;

                //std::cerr<<Xt[0]<<" "<<Xt[1]<<" "<<U[0]<<" "<<Xtp1[0]<<" "<<Xtp1[1]<<" "<<pos_it->second<<" "<<spd_it->second<<std::endl;
                }

            }
        }
        
    }
}

//vector<double> Lipm2Description::learnParams(Value& Xt, Value& U, unsigned int nb_sample_per_dim)
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

Lipm2Description::~Lipm2Description()
{
        //dtor
}


//???
void Lipm2OneGaussDescription::getTransDistrib(Distributions::Distrib& P, Value& Xt, Value& U)
{
    cXt.continuize(Xt);
    cU.continuize(U);
    sim.setXU(cXt,cU);
    cXtp1=sim.simulate(DT,NB_STEP_SIM);

    // Distributions::Distrib1D p_pos, p_spd;
    Distributions::Distrib1D p_pos, p_spd, f_pos, f_spd;
    
    Distributions::dirac(p_pos, state, 0, cXt[0]+DT*cXtp1[1]);
    Distributions::gaussian(p_spd, state, 1, cXtp1[1], SIGMA_S, NB_SIG);
    // Distributions::Distrib1D::iterator pos_it, spd_it;

    Distributions::Distrib1D::iterator pos_it, spd_it, fpos_it, fspd_it;

    Value Xtp1(state);
    for(pos_it=p_pos.begin(); pos_it!=p_pos.end(); ++pos_it)
    {
        for(spd_it=p_spd.begin(); spd_it!=p_spd.end(); ++spd_it)
        {
            Xtp1[0]=pos_it->first;
            Xtp1[1]=spd_it->first;
            P[Xtp1]=pos_it->second * spd_it->second;
        }
    }
}

void Lipm2Simulator::simulate(ContValue& X, ContValue& U, double dt)
{

    double xt=0.0;
    double dxt=0.0;

    double fxt=0.0;
    double dfxt=0.0;





    if(U[1]==1.0) //changement d'appui
    {


        
        // if(U[0]==X[2] && X[3]==X[1]) //pied arriv� (m�me vitesse relative)
        if(U[0]==X[2] && X[3]==X[1]) //pied arriv� (m�me vitesse relative) TODO marge?
        {
        
            // xt=0.5*(X[0]+TAU*X[1])*exp(dt/TAU)+0.5*(X[0]-TAU*X[1])*exp(-dt/TAU);
            // xt=U[0];
            // dxt=0.5/TAU*(X[0]+TAU*X[1])*exp(dt/TAU)-0.5/TAU*(X[0]-TAU*X[1])*exp(-dt/TAU);

            fxt=-X[0]; //switch ancien appui
            dfxt=X[1];
            
            xt=0.5*(U[0]+TAU*X[1])*exp(dt/TAU)+0.5*(U[0]-TAU*X[1])*exp(-dt/TAU);
            dxt=0.5/TAU*(U[0]+TAU*X[1])*exp(dt/TAU)-0.5/TAU*(U[0]-TAU*X[1])*exp(-dt/TAU);
        }
        else 
        {
                //MAJ jambe en vol
            
            fxt=X[2]+X[3]*dt;
            dfxt=X[3]+FOOT_ACC*dt;
            if(dfxt>MAX_FOOT_VEL)
                dfxt=MAX_FOOT_VEL;
        
            xt=0.5*(X[0]+TAU*X[1])*exp(dt/TAU)+0.5*(X[0]-TAU*X[1])*exp(-dt/TAU);
            dxt=0.5/TAU*(X[0]+TAU*X[1])*exp(dt/TAU)-0.5/TAU*(X[0]-TAU*X[1])*exp(-dt/TAU);
            

        }
        
        
        
        
        
    }
    else
    {

        fxt=X[0];
        dfxt=X[1];
        
        xt=0.5*(X[0]+TAU*X[1])*exp(dt/TAU)+0.5*(X[0]-TAU*X[1])*exp(-dt/TAU);
        dxt=0.5/TAU*(X[0]+TAU*X[1])*exp(dt/TAU)-0.5/TAU*(X[0]-TAU*X[1])*exp(-dt/TAU);
            
            
    }
        
    
    
    
    X[0]=xt;
    X[1]=dxt;
    X[2]=fxt;
    X[3]=dfxt;


        //for (unsigned int i=0; i<X.getVar().nbDim(); i++)
        //X.circularize(i);
}


void Lipm2Manager::buildVariables()
{
    vector<unsigned int> stateCards(4);
    vector<double> stateMins(4);
    vector<double> stateMaxs(4);
    vector<bool> stateCircular(4);
    stateCards[0]=XCARD; stateMins[0]=XMIN; stateMaxs[0]=XMAX; stateCircular[0]=false;
    stateCards[1]=VCARD; stateMins[1]=VMIN; stateMaxs[1]=VMAX; stateCircular[1]=false;
    stateCards[2]=FXCARD; stateMins[2]=FXMIN; stateMaxs[2]=FXMAX; stateCircular[2]=false;
    stateCards[3]=FVCARD; stateMins[3]=FVMIN; stateMaxs[3]=FVMAX; stateCircular[3]=false;


    
    state = new Variable(stateMins, stateMaxs, stateCards, stateCircular);

    vector<unsigned int> actionCards(2);
    vector<double> actionMins(2);
    vector<double> actionMaxs(2);
    vector<bool> actionCircular(2);
    actionCards[0]=UCARD; actionMins[0]=UMIN; actionMaxs[0]=UMAX; actionCircular[0]=false;
    actionCards[1]=2; actionMins[1]=0; actionMaxs[1]=1; actionCircular[1]=false;

    action= new Variable(actionMins, actionMaxs, actionCards, actionCircular);


}

void Lipm2Manager::buildSimulator()
{
    sim = new Lipm2Simulator(*state, *action, DT, XMIN, XMAX);
}

void Lipm2Manager::buildModel()
{
    model = new Lipm2Description(*state, *action, *sim);
}

void Lipm2Manager::buildOptimizer(bool drawnPolicy, double initBeta)
{
    optimizer = new OptimizerQuasiDist(*model, drawnPolicy, initBeta);
}

void Lipm2Manager::buildDistribLearner()
{
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 0, (*state)^(*action), NB_SIG, *sim, SIGMA_P, LEARN_INIT_WEIGHT));
    distrib1DLearners.push_back(new CndGaussian1DLearner(*state, 1, (*state)^(*action), NB_SIG, *sim, SIGMA_S, LEARN_INIT_WEIGHT));
}

void Lipm2Manager::buildModelLearner()
{
    learner = new ModelLearner(*optimizer, distrib1DLearners, *state, (*state)^(*action), *sim);
}

void Lipm2Manager::buildController()
{
    controller = new Controller(*optimizer, *learner);
}

void Lipm2Manager::buildObjective()
{
    (*obj_c)[0]=OBJ_P;
    (*obj_c)[1]=OBJ_S;
    (*obj_c)[2]=OBJ_FP;
    (*obj_c)[3]=OBJ_FS;
    
}

void Lipm2Manager::run()
{
//    cout<<"Learning ..."<<endl;
//    learner->learnExhaustive(100,100, DT);
//    learner->update(*obj);
//    cout<<"Done."<<endl;
//    optimizer->computePolicy(*obj);
//    optimizer->toPlot(*obj);
//    save("ex_learn.ser");
  //  exit(0);

    unsigned int max_t=100;
    unsigned int max_t_plot=100;
    Value Xt(*state);
    Value Xtp1(*state);
    Value U(*action);
    ContValue U_c(*action);
    ContValue Xt_c(*state);
    ContValue Xtp1_c(*state);
    Variable stateAction((*state)^(*action));

    ContValue XtU_c(stateAction);
    vector<double> plotPos, plotSpeed, plotAction, plotAction2;
    vector<unsigned int> plotPosDiscr, plotSpeedDiscr, plotActionDiscr, plotActionDiscr2;

    vector<double> sigmas;  /* standard deviation during the control simulation below */
    sigmas.push_back(0.0) /*position*/;
        //sigmas.push_back(0.1); /* speed */
    sigmas.push_back(0.0); /* speed */
    
           /***     Folder where to store logs ****/
    string folder(".");
    optimizer->setDir(folder);

//    optimizer->toPlot(*obj);

//    optimizer->nbEdgesModif=0;
//    unsigned int nb_runs=1;

        //Xt_c[0]=-M_PI;
        //Xt_c[1]=0.0;


        /*
    Xt_c[0]=0.05;
    Xt_c[1]=0.0;

    
    Xt.discretize(Xt_c);
        */


//    optimizer->incrBeta(261);
//    cout<<"BETA = "<<optimizer->getBeta()<<endl;
//    optimizer->computePolicy(*obj);

//    for (unsigned int nb_sim=0; nb_sim<2; nb_sim++)
//    {
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


//     for (unsigned int t=0; t<max_t; t++)
//             {
//             //    cout<<"t="<<t<<" ; U="<<U[0]<<" ; Uc="<<U_c[0]<<endl;
//                 U=controller->control(Xt, *obj);
//                 plotPosDiscr.push_back(Xt[0]);plotSpeedDiscr.push_back(Xt[1]);plotActionDiscr.push_back(U[0]);plotActionDiscr2.push_back(U[1]);
//                 U_c.continuize(U);
//        //         cout<<t<<" ; ("<<Xt[0]<<", "<<Xt[1]<<", "<<U[0]<<")  ;  ("<<Xt_c[0]<<", "<<Xt_c[1]<<", "<<U_c[0]<<")"<<endl;
//                 plotPos.push_back(Xt_c[0]); plotSpeed.push_back(Xt_c[1]); plotAction.push_back(U_c[0]);plotAction2.push_back(U_c[1]);



//                // dynamic_cast<Lipm2Simulator*>(sim)->setFriction(0.2);
//           //      dynamic_cast<Lipm2Simulator*>(sim)->setSlope(0.18);



//                 sim->setXU(Xt_c, U_c);


//            /***     Deterministic sim ****/
//                 Xtp1_c=sim->simulate(DT,NB_STEP_SIM);
// //                cout<<"deterministic : "<<Xtp1_c[0]<<", "<<Xtp1_c[1]<<endl;

//            /***     Gaussian noised sim ****/
// //                Xtp1_c=sim->simulateGaussianDraw(DT, sigmas,5.0, NB_STEP_SIM);


//                 Xtp1.discretize(Xtp1_c);

                
//            /***     Adding data to learn (if needed) ****/
// //                XtU_c=Xt_c^U_c;
// //                learner->addData(Xtp1_c, XtU_c);


//                 Xt_c=Xtp1_c;
//                 Xt=Xtp1;




//            }
//
//        }


            

           /***     Sim logs ****/

                    /*
                
        ofstream file((folder+"/control.dat").c_str());
        for (unsigned int t=0; t<max_t_plot; t++)
        {
            file<<plotPos[t]<<" "<<plotSpeed[t]<<" "<<plotAction[t]<<" "<<plotAction2[t]<<" "<<endl;
        }
        file.close();

        ofstream file_discr((folder+"/control_discr.dat").c_str());
        for (unsigned int t=0; t<max_t_plot; t++)
        {
            file_discr<<plotPosDiscr[t]<<" "<<plotSpeedDiscr[t]<<" "<<plotActionDiscr[t]<<" "<<plotActionDiscr2[t]<<" "<<endl;
        }

        file_discr.close();
        */

            
           /***     landscape and policy saving for plot ****/
        optimizer->toPlot(*obj);




// }

           /***     Parameter learning from previously stored data (if needed) ****///        learner->update(*obj);
//        learner->update(*obj);
//        optimizer->computePolicy(*obj);
//        optimizer->toPlot(*obj);

}

Lipm2Manager::~Lipm2Manager()
{
    //dtor
}
