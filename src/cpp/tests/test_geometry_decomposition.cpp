#include<vector>
#include<iostream>
#include<fstream>
#include<math.h>
#include<geometry_decomposition.h>

#define BOOST_TEST_MODULE test_geometry_decomposition
#include <boost/test/included/unit_test.hpp>

typedef double floatType;
typedef std::vector< floatType > vectorType;
typedef std::vector< vectorType > matrixType;

//bool fuzzy_equals(double a, double b, double tolr=1e-6, double tola=1e-6){
//    /*!
//    Compare two doubles to determine if they are equal.
//    */
//
//    floatType tol = fmin(tolr*fabs(a) + tola, tolr*fabs(b) + tola);
//    return fabs(a-b)<tol;
//}
//
//bool fuzzy_equals(vectorType a, vectorType b, double tolr=1e-6, double tola=1e-6){
//    /*!
//    Compare two vectors to determine if they are equal
//    */
//
//    if (a.size() != b.size()){
//        std::cout << "Error: vectors must have the same size.\n";
//        assert(1==0);
//    }
//
//    for (unsigned int i=0; i<a.size(); i++){
//        if (!fuzzy_equals(a[i], b[i], tolr, tola)){
//            return false;
//        }
//    }
//    return true;
//}
//
//bool fuzzy_equals(matrixType A, matrixType B, double tolr=1e-6, double tola=1e-6){
//    /*!
//    Compare two matrices to determine if they are equal
//    */
//
//    if (A.size() != B.size()){
//        std::cout << "Error: matrices must have the same size.\n";
//        assert(1==0);
//    }
//
//    for (unsigned int i=0; i<A.size(); i++){
//        if (!fuzzy_equals(A[i], B[i], tolr, tola)){
//            return false;
//        }
//    }
//    return true;
//}
//
//template<typename T>
//bool equals(const T &a, const T &b){
//    /*!
//     * Compare two values for equality
//     */
//     return a == b;
//}
//
//template<typename T>
//bool equals(const std::vector< T > &a, const std::vector< T > &b){
//    /*!
//     * Compare two vectors for equality
//     */
//    unsigned int size = a.size();
//    if (size != b.size()){
//        return false;
//    }
//    for (unsigned int i=0; i<size; i++){
//        if (!equals(a[i], b[i])){
//            return false;
//        }
//    }
//    return true;
//}
//
//template<typename T>
//bool equals(const std::vector< std::vector< T > > &a, const std::vector< std::vector< T > > &b){
//    /*!
//     * Compare two matrices for equality
//     */
//
//    unsigned int size = a.size();
//    if (size != b.size()){
//        return false;
//    }
//    for (unsigned int i=0; i<size; i++){
//        if (!equals(a[i], b[i])){
//            return false;
//        }
//    }
//    return true;
//}

void print(vectorType a){
    /*!
    Print the vector to the terminal
    */

    for (unsigned int i=0; i<a.size(); i++){
        std::cout << a[i] << " ";
    }
    std::cout << "\n";
}

void print(matrixType A){
    /*!
    Print the matrix to the terminal
    */

    for (unsigned int i=0; i<A.size(); i++){
        print(A[i]);
    }
}

BOOST_AUTO_TEST_CASE( testGetTets ){
    /*!
     * Test the creation of a collection of tetrahedra that describe 
     * the volume associated with ordered points on a plane that 
     * represent the boundary of a convex polyhedra and some 
     * center point.
     * 
     */

    matrixType nodes = {{-1, -1, 1},
                        { 1, -1, 1},
                        { 1,  1, 1},
                        {-1,  1, 1}};

    vectorType centroid = {0, 0, 0};

    std::vector< matrixType > tets = gDecomp::getTets(centroid, nodes);

    unsigned int i=0, j=1;
    for (auto tet = tets.begin(); tet!=tets.end(); tet++, i++, j++){
        BOOST_CHECK( vectorTools::fuzzyEquals( ( *tet )[ 0 ], centroid ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( ( *tet )[ 1 ], { 0, 0, 1 } ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( ( *tet )[ 2 ], nodes[ i ] ) );
        if (j >= nodes.size()){
            j = 0;
        }
        BOOST_CHECK( vectorTools::fuzzyEquals( ( *tet )[ 3 ], nodes[ j ] ) );
    }

}

BOOST_AUTO_TEST_CASE( testGetTetVolume ){
    /*!
     * Compute the volume of a tetrahedron
     * 
     */

    matrixType tet = {{ 0, 0, 0},
                      { 1, 0, 0},
                      { 0, 1, 0},
                      { 0, 0, 1}};

    double volume = gDecomp::getTetVolume(tet);
    
    BOOST_CHECK( vectorTools::fuzzyEquals( volume, 1./6 ) );

}

BOOST_AUTO_TEST_CASE( testGetUnitToTetMap ){
    /*!
     * Test the computation of the map between the unit tetrahedra and an 
     * arbitrary tetrahedron.
     * 
     */

    matrixType localNodes = {{ 0, 0, 0},
                             { 1, 0, 0},
                             { 0, 1, 0},
                             { 0, 0, 1}};

    matrixType nodes = {{  1,  4, 2},
                        {  6,  4, 1},
                        { 10,  3, 5},
                        {  3, -1, 4}};

    matrixType A;
    vectorType d;
    gDecomp::getUnitToTetMap(nodes, A, d);

    unsigned int i=0;
    for (auto ln=localNodes.begin(); ln!=localNodes.end(); ln++, i++){
        BOOST_CHECK( vectorTools::fuzzyEquals( nodes[ i ], vectorTools::dot( A, *ln ) + d ) );
    }

}

BOOST_AUTO_TEST_CASE( testGetTetQuadrature ){
    /*!
     * Test the quadrature points for the tetrahedron.
     * 
     * 
     */
    
    matrixType points;
    vectorType weights;

    for (unsigned int order=0; order<4; order++){
        gDecomp::getTetQuadrature(order, points, weights);
    }

}

BOOST_AUTO_TEST_CASE( testFindPointsOnFace ){
    /*!
     * Test the utility which detects if points are on a 
     * surface or not.
     * 
     */

    vectorType normal = {1/std::sqrt(3), 1/std::sqrt(3), 1/std::sqrt(3)};
    vectorType point = {1./3, 1./3, 1./3};

    matrixType points = {{ 1.0, 0, 0},
                         { 0.0, 1, 0},
                         { 0.0, 0, 1},
                         {-1.0, 0, 1},
                         { 1.1, 0, 0}};

    std::vector< unsigned int > answers = {0, 1, 2};

    std::vector< unsigned int > r;
    gDecomp::findPointsOnFace(normal, point, points, r);

    BOOST_CHECK( vectorTools::equals( r, answers ) );
}

BOOST_AUTO_TEST_CASE( testOrderPlanarPoints ){
    /*!
     * Test the utility which returns the indices which 
     * order the incoming points CCW
     * 
     */

    matrixType points = {{ 1.0,  0.0, 0},
                         { 0.0, -0.2, 0},
                         { 0.0,  1.0, 0},
                         {-1.0, -1.0, 0}};

    std::vector< unsigned int > idx;
    gDecomp::orderPlanarPoints(points, idx);
    BOOST_CHECK( vectorTools::equals( idx, { 2, 0, 1, 3 } ) );

}

BOOST_AUTO_TEST_CASE( testGetFacePoints ){
    /*!
     * Test the utility which returns the indices of 
     * the points located on each face.
     * 
     */

    matrixType points = { {-0.000000000, 1.000000000, -0.000000000},
                          {-0.000000000, 0.361803399, -0.000000000},
                          {-0.000000000, 0.500000000, 0.500000000},
                          {-0.000000000, 0.361803399, 0.361803399},
                          {0.500000000, 0.500000000, -0.000000000},
                          {0.361803399, 0.361803399, -0.000000000},
                          {0.361803399, 0.361803399, 0.276393202},
                          {0.276393202, 0.361803399, 0.361803399}};

    std::vector< gDecomp::faceType > faces = {std::pair<vectorType, vectorType>({-1.000000000, 0.000000000, 0.000000000},
                                                                                {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.000000000, -1.000000000, 0.000000000},
                                                                                {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.000000000, 0.000000000, -1.000000000},
                                                                                {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.577350269, 0.577350269, 0.577350269},
                                                                                {1.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.707106781, -0.707106781, 0.000000000},
                                                                                {0.361803399, 0.361803399, 0.138196601}),
                                              std::pair<vectorType, vectorType>({0.000000000, -1.000000000, 0.000000000},
                                                                                {0.138196601, 0.361803399, 0.138196601}),
                                              std::pair<vectorType, vectorType>({0.000000000, -0.707106781, 0.707106781},
                                                                                {0.138196601, 0.361803399, 0.361803399})};

    std::vector< std::vector< unsigned int > > answer = {{  3, 2, 0, 1},
                                                         {},
                                                         {  5, 4, 0, 1},
                                                         {  6, 4, 0, 2, 7},
                                                         {  4, 5, 6},
                                                         {  6, 5, 1, 3, 7},
                                                         {  2, 3, 7}};

    std::vector< std::vector< unsigned int > > indexFaces;
    gDecomp::getFacePoints(faces, points, indexFaces);

    BOOST_CHECK( vectorTools::equals( answer, indexFaces ) );

}

BOOST_AUTO_TEST_CASE( testVolumeToTets ){
    /*!
     * Test the utility to deconstruct a volume into 
     * tetrahedra.
     * 
     */

    matrixType hexPoints = {{-1, -1, -1},
                            { 1, -1, -1},
                            { 1,  1, -1},
                            {-1,  1, -1},
                            {-1, -1,  1},
                            { 1, -1,  1},
                            { 1,  1,  1},
                            {-1,  1,  1}};

    std::vector< gDecomp::faceType > hexFaces = {std::pair<vectorType, vectorType>({ 1, 0, 0}, { 1, 0, 0}),
                                                 std::pair<vectorType, vectorType>({-1, 0, 0}, {-1, 0, 0}),
                                                 std::pair<vectorType, vectorType>({ 0, 1, 0}, { 0, 1, 0}),
                                                 std::pair<vectorType, vectorType>({ 0,-1, 0}, { 0,-1, 0}),
                                                 std::pair<vectorType, vectorType>({ 0, 0, 1}, { 0, 0, 1}),
                                                 std::pair<vectorType, vectorType>({ 0, 0,-1}, { 0, 0,-1})};

    std::vector< matrixType > hexTets;
    gDecomp::volumeToTets(hexFaces, hexPoints, hexTets);

    floatType hexVolume = 0;
    for (auto tet=hexTets.begin(); tet!=hexTets.end(); tet++){
        hexVolume += gDecomp::getTetVolume(*tet);
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( hexVolume, 8. ) );

}

BOOST_AUTO_TEST_CASE( testFindMidpoints ){
    /*!
     * Test the computation of the midpoints between a 
     * point and a collection of points.
     * 
     */

    vectorType p = {1, 2, 3};
    matrixType points = {{3, 4, 5},
                         {6, 1, 4},
                         {1, 2, 3},
                         {5, 1,-1}};

    matrixType midpointsAnswer;
    gDecomp::findMidpoints(p, points, midpointsAnswer);

    matrixType midpointsSolution = {{2.0, 3.0, 4.0},
                                    {3.5, 1.5, 3.5},
                                    {3.0, 1.5, 1.0}};

    BOOST_CHECK( vectorTools::fuzzyEquals( midpointsAnswer, midpointsSolution ) );

}

BOOST_AUTO_TEST_CASE( testFindPointOfIntersection ){
    /*!
     * Test the computation of the point of intersection of three planes
     * 
     */
    
    std::vector< gDecomp::faceType > planes = {std::pair< vectorType, vectorType >({1, 0, 0}, {1.0, 0.5, 0.5}),
                                               std::pair< vectorType, vectorType >({0, 1, 0}, {0.5, 1.0, 0.5}),
                                               std::pair< vectorType, vectorType >({0, 0, 1}, {0.5, 0.5, 1.0})};

    vectorType pointAnswer;
    bool solveFlag;
    gDecomp::findPointOfIntersection(planes, pointAnswer, solveFlag);

    BOOST_CHECK( vectorTools::fuzzyEquals( pointAnswer, { 1, 1, 1 } ) );

    BOOST_CHECK( solveFlag );

    planes = {std::pair< vectorType, vectorType >({1, 0, 0}, {1.0, 0.5, 0.5}),
              std::pair< vectorType, vectorType >({0, 1, 0}, {0.5, 1.0, 0.5}),
              std::pair< vectorType, vectorType >({0, 1, 0}, {0.5, 0.5, 1.0})};

    gDecomp::findPointOfIntersection(planes, pointAnswer, solveFlag);

    BOOST_CHECK( !solveFlag );

}

BOOST_AUTO_TEST_CASE( testFindAllPointsOfIntersection ){
    /*!
     * Test for the utility which finds all of the points of 
     * intersection of a set of planes.
     * 
     */

    std::vector< gDecomp::faceType > hexFaces = {std::pair<vectorType, vectorType>({ 1, 0, 0}, { 1, 0, 0}),
                                                 std::pair<vectorType, vectorType>({-1, 0, 0}, {-1, 0, 0}),
                                                 std::pair<vectorType, vectorType>({ 0, 1, 0}, { 0, 1, 0}),
                                                 std::pair<vectorType, vectorType>({ 0,-1, 0}, { 0,-1, 0}),
                                                 std::pair<vectorType, vectorType>({ 0, 0, 1}, { 0, 0, 1}),
                                                 std::pair<vectorType, vectorType>({ 0, 0,-1}, { 0, 0,-1})};

    matrixType intersectionAnswers = {{ 1, 1, 1},
                                      { 1, 1,-1},
                                      { 1,-1, 1},
                                      { 1,-1,-1},
                                      {-1, 1, 1},
                                      {-1, 1,-1},
                                      {-1,-1, 1},
                                      {-1,-1,-1}};

    matrixType intersectionPoints;
    gDecomp::findAllPointsOfIntersection(hexFaces, intersectionPoints);

    BOOST_CHECK( vectorTools::fuzzyEquals( intersectionPoints, intersectionAnswers ) );

    std::vector< gDecomp::faceType > faces = {std::pair<vectorType, vectorType>({-1.000000000, 0.000000000, 0.000000000}, {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.000000000, -1.000000000, 0.000000000}, {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.000000000, 0.000000000, -1.000000000}, {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.577350269, 0.577350269, 0.577350269}, {1.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.707106781, -0.707106781, 0.000000000}, {0.361803399, 0.361803399, 0.138196601}),
                                              std::pair<vectorType, vectorType>({0.000000000, -1.000000000, 0.000000000}, {0.138196601, 0.361803399, 0.138196601}),
                                              std::pair<vectorType, vectorType>({0.000000000, -0.707106781, 0.707106781}, {0.138196601, 0.361803399, 0.361803399})};

    intersectionAnswers = { {-0.000000000, -0.000000000, -0.000000000},
                            {-0.000000000, -0.000000000, 1.000000000},
                            {-0.000000000, 1.000000000, -0.000000000},
                            {-0.000000000, 0.361803399, -0.000000000},
                            {-0.000000000, 0.361803399, 0.638196601},
                            {-0.000000000, 0.500000000, 0.500000000},
                            {-0.000000000, 0.361803399, 0.361803399},
                            {1.000000000, -0.000000000, -0.000000000},
                            {0.500000000, 0.500000000, -0.000000000},
                            {0.638196601, 0.361803399, -0.000000000},
                            {0.361803399, 0.361803399, -0.000000000},
                            {0.361803399, 0.361803399, 0.276393202},
                            {0.333333333, 0.333333333, 0.333333333},
                            {0.276393202, 0.361803399, 0.361803399},
                            {0.361803399, 0.361803399, 0.361803399}};

    gDecomp::findAllPointsOfIntersection(faces, intersectionPoints);

    BOOST_CHECK( vectorTools::fuzzyEquals( intersectionPoints, intersectionAnswers ) );

}

BOOST_AUTO_TEST_CASE( testIsDuplicate ){
    /*!
     * Test of the utility for detecting duplicates in 
     * collections of points.
     * 
     */

    vectorType v = {1, 2, 3};
    matrixType m = {{2, 3,  4},
                    {5, 6,  7},
                    {8, 9, 10},
                    {11, 12, 13}};

    BOOST_CHECK( !gDecomp::isDuplicate( v, m ) );

    v = {5, 6, 7};
    BOOST_CHECK( gDecomp::isDuplicate( v, m ) );

}

BOOST_AUTO_TEST_CASE( testDetermineInteriorPoints ){
    /*!
     * Test the identification of the points located inside of the 
     * domain.
     * 
     */

    std::vector< gDecomp::faceType > faces = {std::pair<vectorType, vectorType>({-1.000000000, 0.000000000, 0.000000000}, {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.000000000, -1.000000000, 0.000000000}, {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.000000000, 0.000000000, -1.000000000}, {0.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.577350269, 0.577350269, 0.577350269}, {1.000000000, 0.000000000, 0.000000000}),
                                              std::pair<vectorType, vectorType>({0.707106781, -0.707106781, 0.000000000}, {0.361803399, 0.361803399, 0.138196601}),
                                              std::pair<vectorType, vectorType>({0.000000000, -1.000000000, 0.000000000}, {0.138196601, 0.361803399, 0.138196601}),
                                              std::pair<vectorType, vectorType>({0.000000000, -0.707106781, 0.707106781}, {0.138196601, 0.361803399, 0.361803399})};

    matrixType points = { {-0.000000000, -0.000000000, -0.000000000},
                          {-0.000000000, -0.000000000, 1.000000000},
                          {-0.000000000, 1.000000000, -0.000000000},
                          {-0.000000000, 0.361803399, -0.000000000},
                          {-0.000000000, 0.361803399, 0.638196601},
                          {-0.000000000, 0.500000000, 0.500000000},
                          {-0.000000000, 0.361803399, 0.361803399},
                          {1.000000000, -0.000000000, -0.000000000},
                          {0.500000000, 0.500000000, -0.000000000},
                          {0.638196601, 0.361803399, -0.000000000},
                          {0.361803399, 0.361803399, -0.000000000},
                          {0.361803399, 0.361803399, 0.276393202},
                          {0.333333333, 0.333333333, 0.333333333},
                          {0.276393202, 0.361803399, 0.361803399},
                          {0.361803399, 0.361803399, 0.361803399}};

    std::vector< unsigned int > interiorPointsAnswers = {2, 3, 5, 6, 8, 10, 11, 13};

    std::vector< unsigned int > interiorPoints;
    gDecomp::determineInteriorPoints({0.1381966, 0.5854102, 0.1381966}, points, faces, interiorPoints);

    BOOST_CHECK( vectorTools::equals( interiorPointsAnswers, interiorPoints ) );

}

BOOST_AUTO_TEST_CASE( testMidpointsToFaces ){
    /*!
     * Test the mapping of the calculated midpoints to the 
     * faces.
     * 
     */

    std::vector< gDecomp::faceType > tetFaces = {std::pair<vectorType, vectorType>({-1.000000000,  0.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000, -1.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.577350269,  0.577350269,  0.577350269},
                                                                                   { 1.000000000,  0.000000000,  0.000000000})};

    matrixType points = {{0.585410197, 0.138196601, 0.138196601},
                         {0.138196601, 0.138196601, 0.138196601},
                         {0.138196601, 0.138196601, 0.585410197},
                         {0.138196601, 0.585410197, 0.138196601}};

    std::vector< gDecomp::faceType > answerFaces = {std::pair<vectorType, vectorType>({-1.000000000,  0.000000000,  0.000000000},
                                                                                      { 0.361803399,  0.138196601,  0.138196601}),
                                                    std::pair<vectorType, vectorType>({-0.707106781,  0.000000000,  0.707106781},
                                                                                      { 0.361803399,  0.138196601,  0.361803399}),
                                                    std::pair<vectorType, vectorType>({-0.707106781,  0.707106781,  0.000000000},
                                                                                      { 0.361803399,  0.361803399,  0.138196601}),
                                                    std::pair<vectorType, vectorType>({ 1.000000000,  0.000000000,  0.000000000},
                                                                                      { 0.361803399,  0.138196601,  0.138196601}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000,  1.000000000},
                                                                                      { 0.138196601,  0.138196601,  0.361803399}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000,  1.000000000,  0.000000000},
                                                                                      { 0.138196601,  0.361803399,  0.138196601}),
                                                    std::pair<vectorType, vectorType>({ 0.707106781,  0.000000000, -0.707106781},
                                                                                      { 0.361803399,  0.138196601,  0.361803399}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                      { 0.138196601,  0.138196601,  0.361803399}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000,  0.707106781, -0.707106781},
                                                                                      { 0.138196601,  0.361803399,  0.361803399}),
                                                    std::pair<vectorType, vectorType>({ 0.707106781, -0.707106781,  0.000000000},
                                                                                      { 0.361803399,  0.361803399,  0.138196601}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000, -1.000000000,  0.000000000},
                                                                                      { 0.138196601,  0.361803399,  0.138196601}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000, -0.707106781,  0.707106781},
                                                                                      { 0.138196601,  0.361803399,  0.361803399})};

    matrixType midpoints;
    std::vector< gDecomp::faceType > midpointFaces;

    for (unsigned int gpt=0; gpt<4; gpt++){
        gDecomp::findMidpoints(points[gpt], points, midpoints);
        gDecomp::midpointsToFaces(points[gpt], midpoints, midpointFaces);

        for (unsigned int f=0; f<3; f++){
            BOOST_CHECK( !( !vectorTools::fuzzyEquals( answerFaces[ gpt*3 + f ].first, midpointFaces[ f ].first, 1e-5 ) &&
                            !vectorTools::fuzzyEquals( answerFaces[ gpt*3 + f ].second, midpointFaces[ f ].second, 1e-5 ) ) );
        }

    }

}

BOOST_AUTO_TEST_CASE( testGetVolumeSubdomainAsTets ){
    /*!
     * Test the decomposition of a volume's subdomain into 
     * tetrahedra. This requires only the planes which form 
     * the volume and a collection of points inside each of 
     * the subdomains.
     * 
     */

    std::vector< gDecomp::faceType > tetFaces = {std::pair<vectorType, vectorType>({-1.000000000,  0.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000, -1.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.577350269,  0.577350269,  0.577350269},
                                                                                   { 1.000000000,  0.000000000,  0.000000000})};

    matrixType points = {{0.585410197, 0.138196601, 0.138196601},
                         {0.138196601, 0.138196601, 0.138196601},
                         {0.138196601, 0.138196601, 0.585410197},
                         {0.138196601, 0.585410197, 0.138196601}};

    std::vector< matrixType > subdomainTets;
    
    vectorType subdomainVolumes(points.size(), 0);
    vectorType tetVolumesAnswer = {0.03980327668541683, 0.04725683661041613, 0.03980327668541683, 0.03980327668541683};

    for (unsigned int index=0; index<points.size(); index++){
        gDecomp::getVolumeSubdomainAsTets(index, points, tetFaces, subdomainTets);
        
        for (auto tet=subdomainTets.begin(); tet!=subdomainTets.end(); tet++){
            subdomainVolumes[index] += gDecomp::getTetVolume(*tet);
        }
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( tetVolumesAnswer, subdomainVolumes ) );

    std::vector< gDecomp::faceType > hexFaces = {std::pair<vectorType, vectorType>({ 1, 0, 0}, { 1, 0, 0}),
                                                 std::pair<vectorType, vectorType>({-1, 0, 0}, {-1, 0, 0}),
                                                 std::pair<vectorType, vectorType>({ 0, 1, 0}, { 0, 1, 0}),
                                                 std::pair<vectorType, vectorType>({ 0,-1, 0}, { 0,-1, 0}),
                                                 std::pair<vectorType, vectorType>({ 0, 0, 1}, { 0, 0, 1}),
                                                 std::pair<vectorType, vectorType>({ 0, 0,-1}, { 0, 0,-1})};

    matrixType hexPoints = {{-1, -1, -1},
                            { 1, -1, -1},
                            { 1,  1, -1},
                            {-1,  1, -1},
                            {-1, -1,  1},
                            { 1, -1,  1},
                            { 1,  1,  1},
                            {-1,  1,  1}};
    hexPoints /=std::sqrt(3);

    subdomainVolumes = vectorType(hexPoints.size(), 0);
    vectorType hexVolumesAnswer = {1, 1, 1, 1, 1, 1, 1, 1};

    for (unsigned int index=0; index<hexPoints.size(); index++){
        gDecomp::getVolumeSubdomainAsTets(index, hexPoints, hexFaces, subdomainTets);
        
        for (auto tet=subdomainTets.begin(); tet!=subdomainTets.end(); tet++){
            subdomainVolumes[index] += gDecomp::getTetVolume(*tet);
        }
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( hexVolumesAnswer, subdomainVolumes ) );

}

BOOST_AUTO_TEST_CASE( testMapLocalTetPointsToGlobal ){
    /*!
     * Test the mapping of local tet points to global coordinates
     * 
     */

    matrixType tet = {{1.39293837, 0.57227867, 0.45370291},
                      {1.10262954, 1.43893794, 0.84621292},
                      {1.9615284,  1.36965948, 0.9618638 },
                      {0.78423504, 0.68635603, 1.45809941}};

    matrixType testPoints = {{0.0, 0.0, 0.0},
                             {1.0, 0.0, 0.0},
                             {0.0, 1.0, 0.0},
                             {0.0, 0.0, 1.0},
                             {0.25, 0.25, 0.25}};

    matrixType globalSolutions = {{1.39293837, 0.57227867, 0.45370291},
                                  {1.10262954, 1.43893794, 0.84621292},
                                  {1.9615284,  1.36965948, 0.9618638 },
                                  {0.78423504, 0.68635603, 1.45809941},
                                  {1.31033284, 1.01680803, 0.92996976}};

    matrixType globalAnswers;
    gDecomp::mapLocalTetPointsToGlobal(tet, testPoints, globalAnswers);

    BOOST_CHECK( vectorTools::fuzzyEquals( globalAnswers, globalSolutions ) );

}

BOOST_AUTO_TEST_CASE( testTetIO ){
    /*!
     * Test the ability to write tetrahedra to a file
     * and read tetrahedra from a file.
     * 
     */

    std::vector< matrixType > tets = {{{ 1,  2,  3},
                                      { 4,  5,  6},
                                      { 7,  8,  9},
                                      {10, 11, 12}},
                                     {{13, 14, 15},
                                      {16, 17, 18},
                                      {19, 20, 21},
                                      {22, 23, 24}},
                                     {{25, 26, 27},
                                      {28, 29, 30},
                                      {31, 32, 33},
                                      {34, 35, 36}}};

    std::string fileName("test.tets");
    gDecomp::writeTetsToFile(fileName, tets);
    std::vector< matrixType > readTets;
    BOOST_CHECK( gDecomp::readTetsFromFile( fileName, readTets ) == 0 );

    unsigned int n=0;
    for (auto tet=readTets.begin(); tet!=readTets.end(); tet++, n++){
        BOOST_CHECK( vectorTools::fuzzyEquals( *tet, tets[ n ] ) );
    }

}

BOOST_AUTO_TEST_CASE( testRemoveDuplicateFaces ){
    /*!
     * Test the utility for removing faces which are identical.
     * 
     */

    std::vector< gDecomp::faceType > tetFaces = {std::pair<vectorType, vectorType>({-1.000000000,  0.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000, -1.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.577350269,  0.577350269,  0.577350269},
                                                                                   { 1.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                   { 1.000000000,  1.000000000,  0.000000000})};

    std::vector< gDecomp::faceType > uniqueFaces = {std::pair<vectorType, vectorType>({-1.000000000,  0.000000000,  0.000000000},
                                                                                      { 0.000000000,  0.000000000,  0.000000000}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000, -1.000000000,  0.000000000},
                                                                                      { 0.000000000,  0.000000000,  0.000000000}),
                                                    std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                      { 0.000000000,  0.000000000,  0.000000000}),
                                                    std::pair<vectorType, vectorType>({ 0.577350269,  0.577350269,  0.577350269},
                                                                                      { 1.000000000,  0.000000000,  0.000000000})};
   

    gDecomp::removeDuplicateFaces(tetFaces);

    BOOST_CHECK( tetFaces.size( ) == 4 );

    unsigned int i=0;
    for (auto face=tetFaces.begin(); face!=tetFaces.end(); face++, i++){
        BOOST_CHECK( vectorTools::fuzzyEquals( face->first, uniqueFaces[ i ].first ) );
        BOOST_CHECK( vectorTools::fuzzyEquals( face->second, uniqueFaces[ i ].second ) );
    }

}

BOOST_AUTO_TEST_CASE( testPrint ){
    /*!
     * Test to make sure the print test command works.
     * 
     */

    std::vector< gDecomp::faceType > tetFaces = {std::pair<vectorType, vectorType>({-1.000000000,  0.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000, -1.000000000,  0.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.000000000,  0.000000000, -1.000000000},
                                                                                   { 0.000000000,  0.000000000,  0.000000000}),
                                                 std::pair<vectorType, vectorType>({ 0.577350269,  0.577350269,  0.577350269},
                                                                                   { 1.000000000,  0.000000000,  0.000000000})};

    //Redirect std::cout
    std::stringstream buffer;
    std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
    gDecomp::print(tetFaces);
    std::string r = buffer.str();
    //Reset std::cout
    std::cout.rdbuf( old );

    std::string answer("-1 0 0 | 0 0 0\n0 -1 0 | 0 0 0\n0 0 -1 | 0 0 0\n0.57735 0.57735 0.57735 | 1 0 0\n");

    BOOST_CHECK( std::strcmp( r.c_str( ), answer.c_str( ) ) == 0 );

}
