#ifndef WALKINGCONFIG_H_INCLUDED
#define WALKINGCONFIG_H_INCLUDED

#define X_CIRC false
#define XCARD 31
#define XMIN -0.6
#define XMAX 0.6
#define VCARD 31
#define VMIN -1.5
#define VMAX 1.5
#define UCARD 21
#define UMIN -1.0
#define UMAX 1.0

#define SIGMA_P 10.0 //0.1  //0.05 // 0.13
#define SIGMA_S 10.0 //0.2 //0.05 //0.13  //position, speed
#define DT 0.25
#define NB_SIG 3.0

#define OBJ_P 0.0
#define OBJ_S 0.0

#define DRAWN_POLICY false
#define NB_STEP_SIM 100
#define LEARN_INIT_WEIGHT 0
#define INIT_BETA 0.1
#define INCR_BETA 0.2


#endif // WALKINGCONFIG_H_INCLUDED
