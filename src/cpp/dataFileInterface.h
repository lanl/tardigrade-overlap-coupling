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

    //XDMF Cell Node Counts
    const std::map< unsigned int, unsigned int > cellNodeCount =
            {
                {  1,  1 }, //Polyvertex
                {  2,  0 }, //Polyline ( special case )
                {  3,  0 }, //Polygon ( special case )
                {  4,  3 }, //Triangle
                {  5,  4 }, //Quadralateral
                {  6,  4 }, //Tetrahedron
                {  7,  5 }, //Pyramid
                {  8,  6 }, //Wedge
                {  9,  8 }, //Hexahedron
                { 16,  0 }, //Polyhedron ( special case )
            };


    //Typedefs
    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< unsigned int > uIntVector; //!Define a vector of unsigned ints

//    //Define the make_shared function
//    template<typename T, typename... Args>
//    std::shared_ptr<T> make_shared(Args&&... args) {
//        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
//    }

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
         * filename: ( The path to the datafile )
         * mode: ( The mode that the file is read ) 
         */

        public:

            //Constructors
            dataFileBase( );
            dataFileBase( const YAML::Node &config );
            dataFileBase( const YAML::Node &config, errorOut error );

            //The factory for child classes
            std::shared_ptr<dataFileBase> create( );
            std::shared_ptr<dataFileBase> create( const std::string &type );

            //File read interface functions
            virtual errorOut getNumIncrements( unsigned int &numIncrements ); //Required overload
            virtual errorOut getNumNodes( const unsigned int increment, unsigned int &numNodes ); //Required overload
            virtual errorOut readMesh( const unsigned int increment, floatVector &nodalPositions ); //Required overload
            virtual errorOut getDomainNodes( const unsigned int increment, const std::string domainName,
                                             uIntVector &domainNodes ); //Required overload
            virtual errorOut getNumDomainNodes( const unsigned int increment, const std::string domainName,
                                                unsigned int &numDomainNodes ); //Required overload
            virtual errorOut getSetNames( const unsigned int increment, std::vector< std::string > &setNames ); //Required overload
            virtual errorOut getSolutionData( const unsigned int increment, const std::string &dataName, const std::string &dataType,
                                              floatVector &data ); //Required overload
            virtual errorOut getMeshData( const unsigned int increment,
                                          floatVector &nodePositions, uIntVector &connectivity, uIntVector &connectivityCellIndices,
                                          unsigned int &cellCounts ); //Required overload

            errorOut _error;
            std::string _filename;
            std::string _mode;

        protected:

            YAML::Node _config;

            errorOut connectivityToCellIndices( const unsigned int &nCells, const uIntVector &connectivity,
                                                uIntVector &connectivityCellIndices );
    };

    class XDMFDataFile : public dataFileBase{
        /*!
         * The XDMF data file interface definition
         *
         * TODO: Consider moving this to its own header at some point.
         *       At this time, there is only one IO file type so it isn't
         *       a big deal.
         */

        public:
            //Constructors
            XDMFDataFile( );
            XDMFDataFile( const YAML::Node &configuration );

            //Overloads
            errorOut getNumIncrements( unsigned int &numIncrements );
            errorOut getNumNodes( const unsigned int increment, unsigned int &numNodes );
            errorOut readMesh( const unsigned int increment, floatVector &nodalPositions );
            errorOut getDomainNodes( const unsigned int increment, const std::string domainName,
                                     uIntVector &domainNodes );
            errorOut getNumDomainNodes( const unsigned int increment, const std::string domainName,
                                        unsigned int &numDomainNodes );
            errorOut getSetNames( const unsigned int increment, std::vector< std::string > &setNames );
            errorOut getSolutionData( const unsigned int increment, const std::string &dataName, const std::string &dataType,
                                      floatVector &data );
            errorOut getMeshData( const unsigned int increment,
                                  floatVector &nodePositions, uIntVector &connectivity, uIntVector &connectivityCellIndices,
                                  unsigned int &cellCounts );

        private:
            //Interface Attributes
            shared_ptr< XdmfReader > _reader;
            shared_ptr< XdmfWriter > _writer;
            shared_ptr< XdmfDomain > _domain;

            //Functions
            void _initializeReadMode( );

            errorOut getXDMFGridCollection( const unsigned int gridCollectionNum,
                                            shared_ptr< XdmfGridCollection > &_gridHolder );

            errorOut getUnstructuredGrid( const unsigned int increment,
                                          shared_ptr< XdmfUnstructuredGrid > &unstructuredGrid );
    };
}

#endif
