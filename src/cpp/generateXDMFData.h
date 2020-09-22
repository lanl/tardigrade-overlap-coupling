/*=============================================================================
|                              generateXDMFData.h                             |
===============================================================================
| Functions used to generate XDMF data which can either be used to test the   |
| overlapCoupling class or to create filtering domains.                       |
=============================================================================*/

#ifndef GENERATEXDMFDATA_H
#define GENERATEXDMFDATA_H

#include<dataFileInterface.h>

namespace fileGenerator{

    typedef dataFileInterface::floatType floatType;
    typedef dataFileInterface::floatVector floatVector;
    typedef dataFileInterface::floatMatrix floatMatrix;
    typedef dataFileInterface::uIntType uIntType;
    typedef dataFileInterface::uIntVector uIntTypeVector;
    typedef dataFileInterface::uIntMatrix uIntTypeMatrix;
    typedef dataFileInterface::errorNode errorNode;
    typedef dataFileInterface::errorOut errorOut;

    class fileGenerator{
        /*!
         * A class which can generate XDMF files as outlined in YAML configuration files
         */
    
        public:
    
            fileGenerator( );
    
            fileGenerator( const std::string &yamlFilename );

            const std::unique_ptr< errorNode > &getError( );
    
        protected:
    
            YAML::Node _config;

            std::unique_ptr< errorNode > _error;
    
    };

}

#endif
