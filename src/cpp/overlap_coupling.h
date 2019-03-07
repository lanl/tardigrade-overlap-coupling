/*!
===============================================================================
|                              overlap_coupling.h                             |
===============================================================================
| The header file for the overlap coupling library. This library provides the |
| classes, functions, and methods required to compute the required terms for  |
| the micro/meso-scale to macro-scale coupling following the micromorphic     |
| continuum mechanics framework.                                              |
===============================================================================
*/

#ifndef OVERLAP_COUPLING_H
#define OVERLAP_COUPLING_H

#include<iostream>
#include<vector>
#include<math.h>
#include<assert.h>

#include "element_library.h"

namespace overlap{

    class OverlapCoupling{
        /*!
        The primary class for performing the overlap coupling operation. This 
        object will return the integration weights for volume integrals, area  
        weighted normal vectors for surface integrals, and also constructs 
        the shape function matrix which can be used to construct the projectors.

        The overlap coupling currently assumes that a gauss domain (i.e. the 
        volume associated with a gauss point) cannot be partially filled with 
        DNS. If a node is in a gauss domain, that full gauss domain is treated 
        as the container for the Voronoi tessellation. 
        */


    };

}

#endif
