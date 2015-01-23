#ifndef TODOROVCONFIG_H_INCLUDED
#define TODOROVCONFIG_H_INCLUDED
#define X_CIRC true   //

//#define XCARD 51
//#define XCARD 60

#define XMIN -M_PI //M_PI-0.4
#define XMAX M_PI //M_PI+0.4

//#define XMIN 0.0 //M_PI-0.4
//#define XMAX 2.0*M_PI //M_PI+0.4


//#define VCARD 51
//#define VCARD 60


#define XCARD 51
#define VCARD 51

//#define VMIN -1.5
//#define VMAX 1.5

#define VMIN -2.5
#define VMAX 2.5

#define UCARD 21

//#define UCARD 10

#define UMIN -0.25
#define UMAX 0.25
/*
#define UMIN -0.5
#define UMAX 0.5
*/

//#define SIGMA_P 0.13 //5.0 //0.13
//#define SIGMA_S 0.13 //5.0 //0.13  //position, speed

//ok
/* #define SIGMA_P 0.2 */
/* #define SIGMA_S 0.2 */

//learning
#define SIGMA_P 0.2
#define SIGMA_S 0.2


/* #define SIGMA_P 0.1 */
/* #define SIGMA_S 0.1 */


#define DT 0.25
#define NB_SIG 3.0

#define OBJ_P 0.0
#define OBJ_S 0.0

#define DRAWN_POLICY false
#define NB_STEP_SIM 1
#define LEARN_INIT_WEIGHT 0


/* #define INIT_BETA 0.1 */
#define INIT_BETA 0.1
#define INCR_BETA 0.2


#endif // TODOROVCONFIG_H_INCLUDED
