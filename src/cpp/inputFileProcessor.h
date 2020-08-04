/*=============================================================================
|                                inputFileProcessor                           |
===============================================================================
| Process input files to put them in a format that can be read by the overlap |
| coupling toolchain.                                                         |
=============================================================================*/

#ifndef INPUTFILEPROCESSOR_H
#define INPUTFILEPROCESSOR

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<yaml-cpp/yaml.h>
#include<unordered_map>

#include<dataFileInterface.h>

namespace inputFileProcessor{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef unsigned int uIntType; //!Define the unsigned ints
    typedef std::vector< uIntType > uIntVector; //!Define a vector of unsigned ints
    typedef std::vector< std::string > stringVector; //!Define a vector of strings
    typedef std::unordered_map< uIntType, uIntType > DOFMap;

    class dataFileReaderBase;

    class inputFileProcessor {
        /*!
         * The class from the file processor which 
         * enables access to output files for use in 
         * coupling a micromorphic macro-scale to a 
         * subscale. The configuration file is in the
         * following format. Parentheses indicate the
         * defaults.
         *
         * macroscale_definition:
         *     filename: path_to_macroscale_file
         *     mode: ( "read" )
         *     filetype: ("XDMF")
         *
         * microscale_definition:
         *     filename: path_to_microscale_file
         *     mode: ( "read" )
         *     filetype: ( "XDMF" )
         *
         * Use of the function entails the following:
         * reader = inputFileProcessor( YAML_filename, macroDataReader, microDataReader );
         */


        public:

            //Constructors
            inputFileProcessor( );
            inputFileProcessor( const std::string &yamlConfigurationFilename );


            //Destructor
            ~inputFileProcessor( );

            //Functions
            errorOut setConfigurationFilename( const std::string &configurationFilename );
            errorOut getError( ){ return _error; }

            const floatVector* getMicroDensities( );
            const floatVector* getMicroVolumes( );
            const floatVector* getMicroWeights( );
            const floatVector* getMicroDisplacements( );
            const floatVector* getMicroBodyForces( );
            const floatVector* getMicroVelocities( );
            const floatVector* getMicroAccelerations( );
            const floatVector* getMicroStresses( );

            const floatVector* getMacroDisplacements( );
            const floatVector* getMacroDispDOFVector( );
            const floatVector* getMacroVelocities( );
            const floatVector* getMacroAccelerations( );

            const floatVector* getMicroNodeReferencePositions( );
            const floatVector* getMacroNodeReferencePositions( );
            const uIntVector*  getMacroNodeReferenceConnectivity( );
            const uIntVector*  getMacroNodeReferenceConnectivityCellIndices( );

            const stringVector* getFreeMicroDomainNames( );
            const stringVector* getGhostMicroDomainNames( );

            const uIntVector* getFreeMicroSurfaceApproximateSplitCount( );
            const uIntVector* getGhostMicroSurfaceApproximateSplitCount( );

            const stringVector* getFreeMicroSurfaceNames( );
            const stringVector* getGhostMicroSurfaceNames( );
            const stringVector* getNonOverlappedMicroSurfaceNames( );

            const uIntVector* getFreeMacroCellIds( );
            const uIntVector* getGhostMacroCellIds( );

            const uIntVector* getFreeMacroCellMicroDomainCounts( );
            const uIntVector* getGhostMacroCellMicroDomainCounts( );

            const uIntVector* getNonOverlappedMicroNodeIds( );
            const uIntVector* getFreeMicroNodeIds( );
            const uIntVector* getGhostMicroNodeIds( );

            const uIntVector* getFreeMacroNodeIds( );
            const uIntVector* getGhostMacroNodeIds( );

            const stringVector *getFreeMacroDomainNames( );
            const stringVector *getGhostMacroDomainNames( );

            const DOFMap *getMicroGlobalToLocalDOFMap( ); 

            const DOFMap *getMacroGlobalToLocalDOFMap( ); 

            bool computeMicroShapeFunctions( );

            const YAML::Node getCouplingInitialization( );

            YAML::Node getVolumeReconstructionConfig( );

            bool microBodyForceDefined( );
            bool microVelocitiesDefined( );
            bool microAccelerationDefined( );
            bool macroVelocitiesDefined( );
            bool macroAccelerationDefined( );

            //Core initialization routines
            errorOut initializeIncrement( const unsigned int microIncrement, const unsigned int macroIncrement );

            //Attributes
            std::shared_ptr< dataFileInterface::dataFileBase > _macroscale;
            std::shared_ptr< dataFileInterface::dataFileBase > _microscale;

            bool useReconstructedMassCenters( );

        private:
            //Private functions
            void initialize( );

            errorOut initializeFileInterfaces( );
            errorOut initializeCouplingDomains( );
            
            errorOut openConfigurationFile( );
            errorOut openConfigurationFile( const std::string &configurationFilename );
            errorOut setMicroNodeWeights( const unsigned int increment );
            errorOut setSurfaceSets( const unsigned int microIncrement );
            errorOut checkCommonDomainConfiguration( YAML::Node domainConfig,
                                                     uIntVector &macroCellIds,
                                                     uIntVector &macroCellMicroDomainCounts,
                                                     stringVector &macroVolumeNodesets,
                                                     stringVector &microVolumeNodesets,
                                                     uIntVector &microSurfaceDomainCount );

            errorOut checkCommonVolumeToSurfaceMapping( const stringVector &microVolumeNodesets, 
                                                        stringVector &microSurfaceNodesets );

            errorOut checkCouplingInitialization( );

            errorOut checkVolumeReconstructionInitialization( );

            errorOut extractDataFileProperties( std::shared_ptr< dataFileInterface::dataFileBase > &dataFile,
                                                const unsigned int &increment, const stringVector &variableKeys,
                                                const std::string &dataType,
                                                const bool &populateWithNullOnUndefined, const std::string &configurationName,
                                                YAML::Node &configuration, bool &populatedFlag, floatVector &properties );

            errorOut extractMicroNodeDensities( const unsigned int &increment );
            errorOut extractMicroBodyForces( const unsigned int &increment );
            errorOut extractMicroVelocities( const unsigned int &increment );
            errorOut extractMicroAccelerations( const unsigned int &increment );
            errorOut extractMicroNodeVolumes( const unsigned int &increment );
            errorOut extractMicroDisplacements( const unsigned int &increment );
            errorOut extractMicroStresses( const unsigned int &increment );
            errorOut extractReferenceMicroMeshData( const unsigned int &increment );

            errorOut extractMacroDisplacements( const unsigned int &increment );
            errorOut extractMacroDispDOFVector( const unsigned int &increment );
            errorOut extractMacroVelocities( const unsigned int &increment );
            errorOut extractMacroAccelerations( const unsigned int &increment );
            errorOut extractReferenceMacroMeshData( const unsigned int &increment );

            errorOut getUniqueNodesInDomains( const unsigned int &increment,
                                              const std::shared_ptr< dataFileInterface::dataFileBase > &dataFile,
                                              const stringVector &domainNames, uIntVector &uniqueIds );

            errorOut setMicroNodeIndexMappings( const unsigned int &increment );
            errorOut setMacroNodeIndexMappings( const unsigned int &increment );

            template< typename T, typename Iter >
            errorOut removeIndicesFromVector( std::vector< T > & v, Iter begin, Iter end );

            //Private Attributes
            bool _increment_initialized = false;
            unsigned int _current_macroIncrement = 0;
            unsigned int _current_microIncrement = 0;
            errorOut _error;
            std::string _configFilename = "";
            YAML::Node _config;
            YAML::Node _volumeReconstructionConfig;

            floatVector _microDomainWeights;
            floatVector _microDensities;
            floatVector _microVolumes;
            floatVector _microDisplacements;
            floatVector _microNodeReferencePositions;
            uIntVector  _microNodeReferenceConnectivity;
            uIntVector  _microNodeReferenceConnectivityCellIndices;
            unsigned int _microCellCounts;

            floatVector _microBodyForces;
            floatVector _microVelocities;
            floatVector _microAccelerations;
            floatVector _microStresses;
            bool _microBodyForceFlag = false;
            bool _microVelocityFlag = false;
            bool _microAccelerationFlag = false;

            floatVector _macroDisplacements;
            floatVector _macroDispDOFVector;
            floatVector _macroVelocities;
            floatVector _macroAccelerations;
            floatVector _macroNodeReferencePositions;
            uIntVector  _macroNodeReferenceConnectivity;
            uIntVector  _macroNodeReferenceConnectivityCellIndices;
            unsigned int _macroCellCounts;
            bool _macroVelocityFlag = false;
            bool _macroAccelerationFlag = false;

            stringVector _non_overlapped_micro_surface_sets;
            stringVector _free_micro_surface_sets;
            stringVector _ghost_micro_surface_sets;

            stringVector _free_micro_volume_sets;
            stringVector _ghost_micro_volume_sets;
            
            uIntVector _free_micro_surface_approximate_split_count;
            uIntVector _ghost_micro_surface_approximate_split_count;

            stringVector _free_macro_volume_sets;
            stringVector _ghost_macro_volume_sets;

            uIntVector _free_macro_cell_ids;
            uIntVector _ghost_macro_cell_ids;

            uIntVector _free_macro_cell_micro_domain_counts;
            uIntVector _ghost_macro_cell_micro_domain_counts;

            uIntVector _unique_non_overlapped_micro_nodes;
            uIntVector _unique_free_micro_nodes;
            uIntVector _unique_ghost_micro_nodes;

            uIntVector _unique_free_macro_nodes;
            uIntVector _unique_ghost_macro_nodes;

            DOFMap _global_to_local_macro_node_map;

            DOFMap _global_to_local_micro_node_map;

            bool _computeMicroShapeFunctions = false;

            uIntType _defaultNumberOfMicroDomainSurfaceRegions = 6;

    };

}

#endif
