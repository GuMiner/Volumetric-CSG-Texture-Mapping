#pragma once
/**************************************************************************/
/* File:   optmem.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

namespace netgen
{

/** 
    Optimized Memory allocation classes
*/
class BlockAllocator
{
private:
  ///
  unsigned size, blocks;
  ///
  void * freelist;
  ///
  Array<char*> bablocks;
public:
  ///
  BlockAllocator (unsigned asize, unsigned ablocks = 100);
  ///
  ~BlockAllocator ();
  ///

  void * Alloc ();

  ///
  void Free (void * p)
  {
    *(void**)p = freelist;
    freelist = p;
  }
};
}