/*=============================================================================
|                              generateXDMFData.h                             |
===============================================================================
| Functions used to generate XDMF data which can either be used to test the   |
| overlapCoupling class or to create filtering domains.                       |
=============================================================================*/

#ifndef GENERATEXDMFDATA_H
#define GENERATEXDMFDATA_H

#include<dataFileInterface.h>
#define USE_EIGEN
#include<vector_tools.h>

namespace fileGenerator{

    typedef dataFileInterface::floatType floatType;
    typedef dataFileInterface::floatVector floatVector;
    typedef dataFileInterface::floatMatrix floatMatrix;
    typedef dataFileInterface::uIntType uIntType;
    typedef dataFileInterface::uIntVector uIntVector;
    typedef dataFileInterface::uIntMatrix uIntMatrix;
    typedef dataFileInterface::errorNode errorNode;
    typedef dataFileInterface::errorOut errorOut;
    typedef dataFileInterface::stringVector stringVector;

    class fileGenerator{
        /*!
         * A class which can generate XDMF files as outlined in YAML configuration files
         */
    
        public:
    
            fileGenerator( );
    
            fileGenerator( const std::string &yamlFilename );

            const std::unique_ptr< errorNode > &getError( );

            int build( );

            const uIntType* getCurrentIncrement( );
    
        protected:
    
            YAML::Node _config;

            std::unique_ptr< errorNode > _error;
            std::shared_ptr< dataFileInterface::dataFileBase > _writer;

        private:

            errorOut _initializeIncrement( const YAML::Node &increment );

            errorOut _writeMeshInformation( const YAML::Node &increment );

            errorOut _writeSolutionInformation( const YAML::Node &increment );

            template< class T >
            errorOut _getPropertyFromYAML( const YAML::Node &node, const std::string &property_name, T &property );

            template< class T >
            errorOut _getKeyValuePairsFromYAML( const YAML::Node &node, const std::string &property_name,
                                                stringVector &keys, std::vector< T > &values );

            uIntType _collectionNumber = 0;

            uIntType _currentIncrement;
    
    };

}

#endif
