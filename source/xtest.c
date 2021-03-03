

/***********************************************************/
/** @file  xtest.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  A test routine to check a part of the program 
 * for any special test.  This is not intended to be 
 * used for anything permanent.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      Conduct a diagnostic test for some purpose
 *
 *
 * @details
 *  This is a totally uncommented routine for use 
 *  with Python to perform some specific diagnostic 
 *  test
 *
 * ### Notes ###
 *
 * xtest is called based on a commmand line switch
 * -xtest.  Depending on whether an exit command is
 * included here (the normal situation), the code
 * will continue to run after that.
 *
 * The routine is called after all of the inputs
 * of a run have been read in.  In many cases
 * (namely when the -r (restart option is used,
 * a previous wind file.
 *
 * xtest is intended to be temporary.  Developers should
 * not expect to be the same the next time they want to
 * debug something.  As a result, developers may want
 * to save their version of xtest.c outside of the
 * git repository for later use.
 *
 *
 *
 **********************************************************/

int
xtest ()
{
//OLD  struct photon p;
//OLD  double s;
//OLD  int allow_negative;
//OLD  int hit;
//OLD  // -333708682.51580548, -3464297957.5598879, 916107434.79302168}, lmn = {-0.3821131382229368, -0.85766764329168899, 0.28692653993169387
//OLD// -3.337e+08 -3.464e+09  9.161e+08
//OLD// -3.821e-01 -8.577e-01  2.869e-01
//OLD  p.x[0] = -333708682.51580548;
//OLD  p.x[1] = -3464297957.5598879;
//OLD  p.x[2] = 916107434.79302168;
//OLD  p.lmn[0] = -0.3821131382229368;
//OLD  p.lmn[1] = -0.85766764329168899;
//OLD  p.lmn[2] = 0.28692653993169387;

//  ([-3236723814.5135012, -5445429831.4690027, -946014556.8735044],
//  [-0.83113545374973818, -0.44527507240102326, -0.9559565373701183])

//OLD  p.x[0] = -3236723814.5135012;
//OLD  p.x[1] = -5445429831.4690027;
//OLD  p.x[2] = -946014556.8735044;
//OLD  p.lmn[0] = -0.83113545374973818;
//OLD  p.lmn[1] = -0.44527507240102326;
//OLD  p.lmn[2] = -0.9559565373701183;

//OLD  allow_negative = 1;

//OLD  s = ds_to_disk (&p, allow_negative, &hit);

//OLD  Log ("s %e\n", s);


  double x[3];
  x[0] = 4.02e+10;
  x[1] = 0;
  x[2] = 1.24e+09;
  double v[3];


  Log ("Start  %.5e %.5e %.5e\n", x[0], x[1], x[2]);
  kn_velocity (0, x, v);

  Log ("Result %.5e %.5e %.5e\n", v[0], v[1], v[2]);



  exit (0);

}
