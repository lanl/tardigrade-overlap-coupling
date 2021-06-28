#include<string>

#include<error_tools.h>
#include<overlapCoupling.h>

typedef overlapCoupling::errorNode errorNode;
typedef overlapCoupling::errorOut errorOut;

int main( int argc, char *argv[] ){

    if ( argc == 1 ){

        errorTools::Node( __func__, "No input file defined. Provide the YAML configuration file.").print( );
	return 1;

    }
    else if ( argc > 2){

        errorTools::Node( __func__, "Too many files defined" );

    }

    // Set the filename
    std::string filename( argv[ 1 ] );

    std::cout << "Constructing overlap coupling object\n";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
	return 1;
    }

    std::cout << "Initializing the overlap coupling object\n";
    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
	return 1;
    }

    unsigned int i = 0;

    std::cout << "Beginning to process increments\n";
    while ( !error ){

	std::cout << "  Processing increment " << i << "\n";
        error = oc.processIncrement( i, 0 );

	if ( error ){

	    std::cerr << "Error in increment i: " << i << "\n";
            error->print( );
	    return 1;

	}

	i++;

    }

    return 0;

}
