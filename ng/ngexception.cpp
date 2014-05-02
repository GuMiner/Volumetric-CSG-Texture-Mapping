/**************************************************************************/
/* File:   ngexception.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 02                                                    */
/**************************************************************************/
#include "stdafx.h"
#include "myadt.hpp"

namespace netgen
{
  NgException :: NgException (const string & s) 
    : what(s)
  { }

  /// append string to description
  void NgException :: Append (const string & s)
  { 
    what += s; 
  }
}
