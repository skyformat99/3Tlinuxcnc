#include "rtapi_math.h"
#include "posemath.h"
#include "7dofkins.h"
#include "kinematics.h"             /* decls for kinematicsForward, etc. */


#include "rtapi.h"		/* RTAPI realtime OS API */
#include "rtapi_app.h"		/* RTAPI realtime module decls */
#include "hal.h"
//#include "motion_debug.h"

struct haldata {
    hal_float_t *A3, *A7, *D2, *D4, *D5, *D7;
    hal_bit_t *INVMOD ;
} *haldata = 0;

#define a3 (*(haldata->A3))
#define a7 (*(haldata->A7))
#define d2 (*(haldata->D2))
#define d4 (*(haldata->D4))
#define d5 (*(haldata->D5))
#define d7 (*(haldata->D7))

static double deltady;
static bool invmod = 0;
int kinematicsForward(const double * joint,
                      EmcPose * world,
                      const KINEMATICS_FORWARD_FLAGS * fflags,
                      KINEMATICS_INVERSE_FLAGS * iflags)
{

   double d1, s2, s3, s4, s5, s6, s7;
   double c2, c3, c4, c5, c6, c7;
   double s34;
   double c34;
   double t1, t2, t3, t4, t5, t6;
   double sumSq, k1, k2, k;
   double px, py, pz;
   PmHomogeneous hom;
   PmPose worldPose;
   PmRpy rpy;

   /* Calculate sin of joints for future use */
   d1 = joint[0];
   s2 = sin(joint[1]*PM_PI/180);
   s3 = sin(joint[2]*PM_PI/180);
   s4 = sin(joint[3]*PM_PI/180);
   s5 = sin(joint[4]*PM_PI/180);
   s6 = sin(joint[5]*PM_PI/180);
   s7 = sin(joint[6]*PM_PI/180);
   /* Calculate cos of joints for future use */
   c2 = cos(joint[1]*PM_PI/180);
   c3 = cos(joint[2]*PM_PI/180);
   c4 = cos(joint[3]*PM_PI/180);
   c5 = cos(joint[4]*PM_PI/180);
   c6 = cos(joint[5]*PM_PI/180);
   c7 = cos(joint[6]*PM_PI/180);
   
   s34 = s3*c4+s4*c3;
   c34 = c3*c4-s3*s4;

   /* Calculate terms to be used in definition of... */
   /* first column of rotation matrix.               */
   t1 = c5*s2 - s5*(c2*c3*c4 - c2*s3*s4);
   t2 = (s2*s5 + c5*(c2*c3*c4 - c2*s3*s4));
   t3 = (c2*c3*s4 + c2*c4*s3);
   t4 = c34*s6 + s34*c5*c6;
   t5 = s2*s3*s4 - c3*c4*s2;
   t6 = c3*s2*s4 + c4*s2*s3;


   /* Define first column of rotation matrix */
   hom.rot.x.x =  s7*(t1) + c7*(c6*t2 - s6*t3);
   hom.rot.x.y = - s7*(c2*c5 - s5*(t5)) - c7*(c6*(c2*s5 + c5*(t5)) + s6*(t6));
   hom.rot.x.z = c7*(t4) - s34*s5*s7;

   /* Define second column of rotation matrix */
   hom.rot.y.x = c7*(t1) - s7*(c6*t2 - s6*t3);
   hom.rot.y.y = s7*(c6*(c2*s5 + c5*(t5)) + s6*(t6)) - c7*(c2*c5 - s5*(t5));
   hom.rot.y.z = - s7*(t4) - s34*c7*s5;


   /* Define third column of rotation matrix */
   hom.rot.z.x = s6*t2 + c6*t3;
   hom.rot.z.y = c6*(t6) - s6*(c2*s5 + c5*(t5));
   hom.rot.z.z = s34*c5*s6 - c34*c6;
   px = d4*s2 + d5*t3 + a3*c2*c3;
   py = d1 + d5*(t6) - d4*c2   + a3*c3*s2;
   pz = d2 - d5*c34 + a3*s3;
	
   /* Define position vector */
//   hom.tran.x = d4*s2 + d7*(s6*t2 + c6*t3) + d5*t3 + a3*c2*c3 + a7*s7*(t1) + a7*c7*(c6*t2 - s6*t3);
//  hom.tran.y = d1 + d5*(t6) - d4*c2 - d7*(s6*(c2*s5 + c5*(t5)) - c6*(t6)) + a3*c3*s2 - a7*s7*(c2*c5 - s5*(t5)) - a7*c7*(c6*(c2*s5 + c5*(t5));
//   hom.tran.z = d2 - d5*c34 + a3*s3 - d7*c34*c6 + a7*c7*(t4) + d7*s34*c5*s6 - a7*s34*s5*s7;
   hom.tran.x = px + a7*hom.rot.x.x + d7*hom.rot.z.x;
   hom.tran.y = py + a7*hom.rot.x.y + d7*hom.rot.z.y;
   hom.tran.z = pz + a7*hom.rot.x.z + d7*hom.rot.z.z;
   // calculate terms for inverse kinematic
   deltady = s2*(d5*s34 + a3*c3) - d4*c2 ;

   /* Calculate terms to be used to...   */
   /* determine flags.                   */
   sumSq = (d1-py)*(d1-py)+ px*px;
  
   /* reset flags */
   *iflags = 0;

   /* Set shoulder-up flag if necessary */
   if (fabs(joint[1]*PM_PI/180 + atan2(d1-py, px) -
       atan2(d4, -sqrt(sumSq-d4*d4))) < FLAG_FUZZ)
   {
     *iflags |= SHOULDER_RIGHT;
   }
   k1 = px*c2 - d1*s2 + py*s2;
   k2 = pz-d2;
   k = (k1*k1+k2*k2-d5*d5-a3*a3)/(2*a3*d5);
   /* Set elbow down flag if necessary */
   if (fabs(joint[3]*PM_PI/180 - atan2(k, sqrt(1-k*k))) < FLAG_FUZZ)
   {
      *iflags |= ELBOW_DOWN;
   }

   /* set singular flag if necessary */
   t1 = hom.rot.z.x*c2*c34+ hom.rot.z.y* s2*c34 + hom.rot.z.z*s34;
   t2 = hom.rot.z.x*s2 - hom.rot.z.y*c2 ;
   if (fabs(t1) < SINGULAR_FUZZ && fabs(t2) < SINGULAR_FUZZ)
   {
      *iflags |= SINGULAR;
   }

   /* if not singular set wrist flip flag if necessary */
   else{
     if (!(fabs(joint[4]*PM_PI/180 - atan2(t2, t1)) < FLAG_FUZZ))
     {
       *iflags |= WRIST_FLIP;
     }
   }

   /* convert hom.rot to world->quat */
   pmHomPoseConvert(&hom, &worldPose);
   pmQuatRpyConvert(&worldPose.rot,&rpy);
   world->tran = worldPose.tran;
   world->a = rpy.r * 180.0/PM_PI;
   world->b = rpy.p * 180.0/PM_PI;
   world->c = rpy.y * 180.0/PM_PI;

   
   /* return 0 and exit */
   return 0;
}


int kinematicsInverse(const EmcPose * world,
                      double * joint,
                      const KINEMATICS_INVERSE_FLAGS * iflags,
                      KINEMATICS_FORWARD_FLAGS * fflags)
{
		invmod = *haldata->INVMOD;

	   PmHomogeneous hom;
	   PmPose worldPose;
	   PmRpy rpy;
	   int singular;

	   double t1, t2;
	   double k1, k2, k;
	   double sumSq;
	   double px, py, pz;

	  double d1;
	   double th3;
	   double th2;
	   double th4;
	   double th5;
	   double th6;
	   double th7;
	   double s3, c3;
	   double s2, c2;
	   double s4, c4;
	   double s5, c5;
	   double s6, c6;
	   double s7, c7, c34, s34;


	   /* reset flags */
	   *fflags = 0;

	   /* convert pose to hom */
	   worldPose.tran = world->tran;
	   rpy.r = world->a*PM_PI/180.0;
	   rpy.p = world->b*PM_PI/180.0;
	   rpy.y = world->c*PM_PI/180.0;
	   pmRpyQuatConvert(&rpy,&worldPose.rot);
	   pmPoseHomConvert(&worldPose, &hom);
	   px = hom.tran.x - a7*hom.rot.x.x - d7*hom.rot.z.x;
	   py = hom.tran.y - a7*hom.rot.x.y - d7*hom.rot.z.y;
	   pz = hom.tran.z - a7*hom.rot.x.z - d7*hom.rot.z.z;
	   /* Joint 1 (2 independent solutions) */
	   if(invmod)
	   d1 = py - deltady;
	   else 
	   d1 = joint[0];
	   // joint 2
	   /* save sum of squares for this and subsequent calcs */
	   sumSq = (d1-py)*(d1-py)+ px*px;

	   /* FIXME-- is use of + sqrt shoulder right or left? */
	   if (*iflags & SHOULDER_RIGHT){
	     th2 = -atan2(d1-py, px) + atan2(d4, -sqrt(sumSq-d4*d4));
	   }
	   else{
	     th2 = -atan2(d1-py, px) + atan2(d4, sqrt(sumSq-d4*d4));
	   }

	   /* save sin, cos for later calcs */
	   s2 = sin(th2);
	   c2 = cos(th2);

	   /* Joint 4 (2 independent solutions) */
	   k1 = px*c2 - d1*s2 + py*s2;
	   k2 = pz-d2;
	   k = (k1*k1+k2*k2-d5*d5-a3*a3)/(2*a3*d5);
	/*   if(k>=1 || k<=-1)
		{	
	//	reportError(_("robot reach limit coordinate"));
		*fflags |= ROBOT_REACH;
		return 0;
		}
*/
	   /* FIXME-- is use of + sqrt elbow up or down? */
	   if (*iflags & ELBOW_DOWN){
	     th4 = atan2(k, sqrt(1-k*k));
	   }
	   else{
	     th4 = atan2(k, -sqrt(1-k*k));
	   }

	   /* compute sin, cos for later calcs */
	   s4 = sin(th4);
	   c4 = cos(th4);

	   /* Joint 3 */

	   c3 = (k1*(d5*s4+a3)-d5*c4*k2)/(d5*c4*d5*c4+(d5*s4+a3)*(d5*s4+a3));
	   s3 = (k2+d5*c3*c4)/(a3+d5*s4) ;

	   th3 = atan2(s3,c3);
	   /* compute sin, cos for later calcs */
	   s3 = sin(th3);
	   c3 = cos(th3);
	//----------------------------------------------------------
	   /* Joint 5 */
		c34 = c3*c4-s3*s4;
		s34 = s3*c4+s4*c3;
	   t1 = hom.rot.z.x*c2*c34+ hom.rot.z.y* s2*c34 + hom.rot.z.z*s34 ;
	   t2 = hom.rot.z.x*s2 - hom.rot.z.y*c2;
	   if (fabs(t1) < SINGULAR_FUZZ && fabs(t2) < SINGULAR_FUZZ){
	     th5 = joint[4]*PM_PI/180;            /* use current value */
	   }
	   else{
		if (fabs(t2) < SINGULAR_FUZZ) t2 = fabs(t2);
	     singular = 0;
	     th5 = atan2(t2, t1);
	   }

	   /* compute sin, cos for later calcs */
	   s5 = sin(th5);
	   c5 = cos(th5);
		th5 = atan2(s5,c5);
	   /* Joint 6 */

	   s6 = hom.rot.z.x*(s2*s5 + c2*c34*c5) + hom.rot.z.y*(s2*c5*c34 - c2*s5) + hom.rot.z.z*s34*c5;
	   c6 = hom.rot.z.x*c2*s34 + hom.rot.z.y*s2*s34 - hom.rot.z.z*c34;
	   th6 = atan2(s6, c6);

	   /* Joint 7 */

	   s7 = hom.rot.x.x*(s2*c5 - c2*c34*s5) - hom.rot.x.y*(s2*c34*s5 + c2*c5) - hom.rot.x.z*s34*s5;
	   c7 = hom.rot.y.x*(s2*c5 - c2*c34*s5) - hom.rot.y.y*(c2*c5 + s2*c34*s5) - hom.rot.y.z*s34*s5;
	   th7 = atan2(s7, c7);

	   /* FIXME-- is wrist flip the normal or offset results? */
	   if (*iflags & WRIST_FLIP){
		if(th5>0)
	     th5 = th5 - PM_PI;
		else
		th5 = th5 + PM_PI;
	     th6 = -th6;
		if(th7>0)
	     th7 = th7 - PM_PI;
		else
		th7 = th7 + PM_PI;
	   }
	   /* copy out */
	   joint[0] = d1;
	   joint[1] = th2*180/PM_PI;
	   joint[2] = th3*180/PM_PI;
	   joint[3] = th4*180/PM_PI;
	   joint[4] = th5*180/PM_PI;
	   joint[5] = th6*180/PM_PI;
		joint[6] = th7*180/PM_PI;
	   return singular == 0 ? 0 : -1;
 	
}
int kinematicsHome(EmcPose * world,
                   double * joint,
                   KINEMATICS_FORWARD_FLAGS * fflags,
                   KINEMATICS_INVERSE_FLAGS * iflags)
{
  /* use joints, set world */
  return kinematicsForward(joint, world, fflags, iflags);
}

KINEMATICS_TYPE kinematicsType()
{
//  return KINEMATICS_FORWARD_ONLY;
  return KINEMATICS_BOTH;
}


EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);

int comp_id;

int rtapi_app_main(void) {
    int res=0;
    
    comp_id = hal_init("7dofkins");
    if (comp_id < 0) return comp_id;
    haldata = hal_malloc(sizeof(struct haldata));
    if (!haldata) goto error;

    if((res = hal_pin_float_new("robot7dof.A3", HAL_IO, &(haldata->A3), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("robot7dof.A7", HAL_IO, &(haldata->A7), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("robot7dof.D2", HAL_IO, &(haldata->D2), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("robot7dof.D4", HAL_IO, &(haldata->D4), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("robot7dof.D5", HAL_IO, &(haldata->D5), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("robot7dof.D7", HAL_IO, &(haldata->D7), comp_id)) < 0) goto error;
    if((res = hal_pin_bit_new("robot7dof.invmod", HAL_IN, &(haldata->INVMOD), comp_id)) < 0) goto error;
    a3 = DEFAULT_A3;
    a7 = DEFAULT_A7;
    d2 = DEFAULT_D2;
    d4 = DEFAULT_D4;
    d5 = DEFAULT_D5;
    d7 = DEFAULT_D7;
    hal_ready(comp_id);
    return 0;
    
error:
    hal_exit(comp_id);
    return res;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }
