/*****************************************************************
* Description: pumakins.h
*   Kinematics for a puma typed robot
*
*   Derived from a work by Fred Proctor
*
* Author: 
* License: GPL Version 2
* System: Linux
*    
* Copyright (c) 2004 All rights reserved.
*
* Last change:
*******************************************************************
* This is the header file to accompany pumakins.c.  
*******************************************************************
*/

#ifndef ROBOT7DOF_H
#define ROBOT7DOF_H

/* the default values for a robot 7dof, these can be changed as HAL parameters */
#define DEFAULT_D2 700
#define DEFAULT_A3 1200
#define DEFAULT_D4 0
#define DEFAULT_D5 430
#define DEFAULT_A7 50
#define DEFAULT_D7 250
#define FLAG_FUZZ 0.000001
#define SINGULAR_FUZZ 0.000001


/* flags for inverse kinematics */
#define SHOULDER_RIGHT 0X01
#define ELBOW_DOWN	0X02
#define WRIST_FLIP	0X04  /* joints at a singularity */
#define SINGULAR	0X08
/* flags for forward kinematics */
#define ROBOT_REACH          0x01  /* pose out of reach */
#endif

