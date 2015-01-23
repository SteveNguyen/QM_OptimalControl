#ifndef LIPMCONFIG_H_INCLUDED
#define LIPMCONFIG_H_INCLUDED
#define X_CIRC false
#define XCARD 21
#define XMIN -0.1
#define XMAX 0.1
#define VCARD 21
#define VMIN -0.1
#define VMAX 0.1

#define FXCARD 5
#define FXMIN -0.1
#define FXMAX 0.1
#define FVCARD 5
#define FVMIN -0.1
#define FVMAX 0.1

#define FOOT_ACC 0.01
#define MAX_FOOT_VEL 0.1

#define UCARD 5
#define UMIN -0.01
#define UMAX 0.01

#define SIGMA_P 0.01 //0.1  //0.05 // 0.13
#define SIGMA_S 0.01 //0.2 //0.05 //0.13  //position, speed
#define DT 0.05
#define NB_SIG 3.0

#define OBJ_P 0.0
#define OBJ_S 0.0
#define OBJ_FP 0.0
#define OBJ_FS 0.0


#define DRAWN_POLICY false
#define NB_STEP_SIM 100
#define LEARN_INIT_WEIGHT 0
#define INIT_BETA 0.1
#define INCR_BETA 0.2

#define TAU 1.0

#endif 
