#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED


//#define _COST_ACC

#define L1 0.2
#define L2 0.2

#define EPSILON 0.00001

#define MINHFOOT 0.03

#define cotan(x) tan(M_PI_2 - x)

#define THETACARD 21
#define DTHETACARD 21
#define HCARD 21
#define DHCARD 21



#define THETAMIN -M_PI/4.0
#define THETAMAX M_PI/4.0

/*
#define DTHETAMIN -M_PI/4.0
#define DTHETAMAX M_PI/4.0
*/

#define DTHETAMIN -0.3
#define DTHETAMAX 0.3

#define HMIN L1
#define HMAX L1+L2

#define DHMIN -0.05
#define DHMAX 0.05




#define DALPHACARD 11
#define DALPHAMIN -0.1
#define DALPHAMAX 0.1

#define DBETACARD 11
#define DBETAMIN -0.1
#define DBETAMAX 0.1

#define SIGMA_P 0.001
#define SIGMA_V 0.001


#define DT 0.25
#define NB_SIG 3.0

#define OBJ_THETA M_PI/8.0
#define OBJ_DTHETA 0.0
#define OBJ_H HMAX-EPSILON //HMAX-0.001
#define OBJ_DH 0.0


#define DRAWN_POLICY false
#define NB_STEP_SIM 1
#define LEARN_INIT_WEIGHT 0
#define INIT_BETA 0.1
#define INCR_BETA 0.2


#endif
