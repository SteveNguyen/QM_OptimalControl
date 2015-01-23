#ifndef PAPERCONFIG_H_INCLUDED
#define PAPERCONFIG_H_INCLUDED
//#define X_CIRC true   //

//#define XCARD 25
#define XCARD 51

//#define XMIN -M_PI //M_PI-0.4
//#define XMAX M_PI //M_PI+0.4


//#define XMIN -4.0 //M_PI-0.4
//#define XMAX 4.0 //M_PI+0.4

#define XMIN -5.0 //M_PI-0.4
#define XMAX 5.0 //M_PI+0.4


//#define XMIN 0.0 //M_PI-0.4
//#define XMAX 2.0*M_PI //M_PI+0.4


//#define YCARD 25
#define YCARD 51


//#define VMIN -1.5
//#define VMAX 1.5

//#define YMIN -4.0
//#define YMAX 4.0

#define YMIN -5.0
#define YMAX 5.0


//#define THETACARD 25
#define THETACARD 51


#define THETAMIN -M_PI
#define THETAMAX M_PI


/*
#define THETAMIN (0.0)
#define THETAMAX (2.0*M_PI)
*/

#define UCARD 11

//#define UCARD 10

#define UMIN -1.0
#define UMAX 1.0

/*
#define UMIN -0.5
#define UMAX 0.5
*/

//#define SIGMA_P 0.13 //5.0 //0.13
//#define SIGMA_S 0.13 //5.0 //0.13  //position, speed

/*
#define SIGMA_P 0.12
#define SIGMA_P 0.12
*/

#define SIGMA_P 0.05
//#define SIGMA_P 0.0


#define DT 0.25
#define NB_SIG 3.0

#define OBJ_PX 0.0
#define OBJ_PY 0.0
#define OBJ_THETA 0.0


#define DRAWN_POLICY true
#define NB_STEP_SIM 1
#define LEARN_INIT_WEIGHT 0
#define INIT_BETA 0.1
#define INCR_BETA 0.2


#endif
