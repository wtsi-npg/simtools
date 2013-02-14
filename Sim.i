/*
 * Import file for SWIG
 *
 * Author: Jennifer Liddle (js10)
 *
 * $Id: Gtc.i 1232 2010-07-22 13:05:56Z js10 $
 *
 */

%module Sim
 %{
 /* Includes the header in the wrapper code */
 #include "Sim.h"
 %}
 
 /* Parse the header file to generate wrappers */
 %include "Sim.h"
