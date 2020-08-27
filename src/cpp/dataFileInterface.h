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
    typedef unsigned int uIntType;
    typedef std::vector< uIntType > uIntVector; //!Define a vector of unsigned ints
    typedef std::vector< uIntVector > uIntMatrix; //!Define a matrix of unsigned ints
    typedef std::vector< std::string > stringVector; //!Define a vector of strings

    //XDMF Cell Node Counts

    const std::map< uIntType, uIntType > cellNodeCount =
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
            virtual errorOut getIncrementTime( const uIntType &increment, floatType &time ); //Required overload
            virtual errorOut getNumIncrements( uIntType &numIncrements ); //Required overload
            virtual errorOut getNumNodes( const uIntType increment, uIntType &numNodes ); //Required overload
            virtual errorOut readMesh( const uIntType increment, floatVector &nodalPositions ); //Required overload
            virtual errorOut getNodeIds( const uIntType increment, const std::string &nodeIdAttributeName, uIntVector &domainNodes );
            virtual errorOut getSubDomainNodes( const uIntType increment, const std::string subDomainName,
                                                uIntVector &domainNodes ); //Required overload
            virtual errorOut getNumSubDomainNodes( const uIntType increment, const std::string subDomainName,
                                                   uIntType &numDomainNodes ); //Required overload
            virtual errorOut getSetNames( const uIntType increment, stringVector &setNames ); //Required overload
            virtual errorOut getSolutionData( const uIntType increment, const std::string &dataName, const std::string &dataType,
                                              floatVector &data ); //Required overload
            virtual errorOut getSolutionVectorDataFromComponents( const uIntType increment,
                                                                  const stringVector &componentNames,
                                                                  const std::string &dataType, floatVector &data ); //Probably doesn't need to be overloaded

            virtual errorOut getMeshData( const uIntType increment,
                                          floatVector &nodePositions, uIntVector &connectivity, uIntVector &connectivityCellIndices,
                                          uIntType &cellCounts ); //Required overload

            virtual errorOut initializeIncrement( const floatType time, const uIntType &reference_increment,
                                                  const uIntType &collectionNumber, uIntType &increment );

            virtual errorOut addRootCollection( const std::string &collectionName, const std::string &collectionDescription,
                                                uIntType &collectionNumber );

            virtual errorOut writeIncrementMeshData( const uIntType increment, const uIntType collectionNumber,
                                                     const uIntVector &nodeIds, const uIntMatrix &nodeSets,
                                                     const stringVector &nodeSetNames, const floatVector &nodePositions,
                                                     const uIntVector &elementIds, const uIntMatrix &elementSets,
                                                     const stringVector &elementSetNames, const uIntVector &connectivity ); //Required overload
            virtual errorOut writeSolutionData( const uIntType increment, const uIntType collectionNumber,
                                                const std::string &dataName, const std::string &dataType,
                                                const floatVector &data );

            errorOut _error;
            std::string _filename;
            std::string _mode;

        protected:

            YAML::Node _config;

            errorOut connectivityToCellIndices( const uIntType &nCells, const uIntVector &connectivity,
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
            XDMFDataFile( YAML::Node configuration );

            //Overloads
            errorOut getIncrementTime( const uIntType &increment, floatType &time );
            errorOut getNumIncrements( uIntType &numIncrements );
            errorOut getNumNodes( const uIntType increment, uIntType &numNodes );
            errorOut readMesh( const uIntType increment, floatVector &nodalPositions );
            errorOut getNodeIds( const uIntType increment, const std::string &nodeIdAttributeName, uIntVector &domainNodes );
            errorOut getSubDomainNodes( const uIntType increment, const std::string subDomainName,
                                        uIntVector &domainNodes );
            errorOut getNumSubDomainNodes( const uIntType increment, const std::string subDomainName,
                                           uIntType &numDomainNodes );
            errorOut getSetNames( const uIntType increment, stringVector &setNames );
            errorOut getSolutionData( const uIntType increment, const std::string &dataName, const std::string &dataType,
                                      floatVector &data );
            errorOut getMeshData( const uIntType increment,
                                  floatVector &nodePositions, uIntVector &connectivity, uIntVector &connectivityCellIndices,
                                  uIntType &cellCounts );

            errorOut initializeIncrement( const floatType time, const uIntType &reference_increment,
                                          const uIntType &collectionNumber, uIntType &increment );
            errorOut addRootCollection( const std::string &collectionName, const std::string &collectionDescription,
                                        uIntType &collectionNumber );
            errorOut writeIncrementMeshData( const uIntType increment, const uIntType collectionNumber, 
                                             const uIntVector &nodeIds, const uIntMatrix &nodeSets,
                                             const stringVector &nodeSetNames, const floatVector &nodePositions,
                                             const uIntVector &elementIds, const uIntMatrix &elementSets,
                                             const stringVector &elementSetNames, const uIntVector &connectivity );

        private:
            //Interface Attributes
            shared_ptr< XdmfReader > _reader;
            shared_ptr< XdmfWriter > _writer;
            shared_ptr< XdmfDomain > _domain;

            uIntVector _increment_reference_grids;

            //Functions
            void _initializeReadMode( );
            void _initializeWriteMode( YAML::Node &configuration );

            errorOut getXDMFGridCollection( const uIntType gridCollectionNum,
                                            shared_ptr< XdmfGridCollection > &_gridHolder );

            errorOut getUnstructuredGrid( const uIntType increment,
                                          shared_ptr< XdmfUnstructuredGrid > &unstructuredGrid );
    };
}

#endif
