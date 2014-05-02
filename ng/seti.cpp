#include "stdafx.h" 
#include "mystdlib.h"
#include "myadt.hpp"

namespace netgen
{
    IndexSet :: IndexSet (int maxind)
    {
        SetMaxIndex (maxind);
    }

    IndexSet :: ~IndexSet ()
    {
        Clear();
    }


    void IndexSet :: SetMaxIndex (int maxind)
    {
        if (maxind > flags.Size())
        {
            flags.SetSize (2 * maxind);
            flags.Clear();
        }
    }

    void IndexSet :: Del (int ind)
    {
        for (int i = 1; i <= set.Size(); i++)
            if (set.Get(i) == ind)
            {
                set.DeleteElement (ind);
                break;
            }
            flags.Clear (ind);
    }

    void IndexSet :: Clear ()
    {
        for (int i = 1; i <= set.Size(); i++)
            flags.Clear (set.Get(i));
        set.SetSize (0);
    }
}
