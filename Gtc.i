/*
 * Import file for SWIG
 *
 * Author: Jennifer Liddle (js10)
 *
 * $Id: Gtc.i 1232 2010-07-22 13:05:56Z js10 $
 *
 */

%module Gtc
 %{
 /* Includes the header in the wrapper code */
 #include "gtc_process.h"
 #include "Gtc.h"
 #include "Manifest.h"
 #include "win2unix.h"
 %}
 
 /* Parse the header file to generate wrappers */
 %include "Gtc.h"
 %include "Manifest.h"
 %include "gtc_process.h"
 %include "win2unix.h"
