/*=============================================================================
|                              dataFileInterface                              |
===============================================================================
| A collection of tools that allow for data-files to be interacted wtih both  |
| for reading and writing of those files.                                     |
===============================================================================
| Currently supported file standards:                                         |
| XDMF                                                                        |
=============================================================================*/

#ifndef DATAFILEINTERFACE_H
#define DATAFILEINTERFACE_H

#include<error_tools.h>
#include<yaml-cpp/yaml.h>
#include<map>
#include<memory>

//XDMF headers

#include "XdmfDomain.hpp"
#include "XdmfInformation.hpp"
#include "XdmfReader.hpp"
#include "XdmfWriter.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfGeometry.hpp"
#include "XdmfGridCollection.hpp"
#include "XdmfGridCollectionType.hpp"
#include "XdmfError.hpp"

namespace dataFileInterface{

    //Typedefs
    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< unsigned int > uIntVector; //!Define a vector of unsigned ints

    //Define the make_unique function
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    //All new readers must be defined both in the registry enum and in the map to allow
    //them to be accessed by strings. They must also be registered in dataFileBase::create

    enum registry { XDMF };
    static std::map< std::string, registry > registryMap = { { "XDMF",  registry::XDMF } };

    class dataFileBase {
        /*!
         * The base class for the datafile
         *
         * Takes in a YAML file which is required to have the following attributes:
         *
         * filename: ( THe 
         */

        public:

            //Constructors
            dataFileBase( );
            dataFileBase( const YAML::Node &config );
            dataFileBase( const YAML::Node &config, errorOut error );

            //The factory for child classes
            std::unique_ptr<dataFileBase> create( );
            std::unique_ptr<dataFileBase> create( const std::string &type );

            //File read interface functions
            virtual errorOut readMesh( const unsigned int increment, floatVector &nodalPositions ); //Required overload

            errorOut _error;
            std::string _filename;
            std::string _mode;

        private:

            YAML::Node _config;
    };

    class XDMFDataFile : public dataFileBase{
        /*!
         * The XDMF data file interface definition
         *
         * TODO: Consider moving this to it's own header at some point
         *       at this time, there is only one IO file type so it isn't
         *       a bit deal.
         */

        public:
            //Constructors
            XDMFDataFile( );
            XDMFDataFile( const YAML::Node &configuration );

            //Overloads
            errorOut readMesh( const unsigned int increment, floatVector &nodalPositions );

        private:
            //Interface Attributes
            shared_ptr< XdmfReader > _reader;
            shared_ptr< XdmfWriter > _writer;
            shared_ptr< XdmfDomain > _domain;

            //Functions
            void _initializeReadMode( );
    };
}

#endif
