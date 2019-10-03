/******************************************************************************
*                          geometry_decomposition.cpp                         *
===============================================================================
* A collection of functions and related utilities intended to help decompose  *
* simple geometric objects into easy to integrate subvolumes.                 *
*******************************************************************************
*/

#include "geometry_decomposition.h"

namespace gDecomp{

    std::vector< matrixType > getTets(const vectorType &p, const matrixType &nodes){
        /*!
         * Get the tetrahedra as lists of points which all use point p
         * as one of the vertices. Nodes should be an ordered list of 
         * points from a flat surface
         * 
         * :param const vectorType &p: Point p
         * :param const matrixType &nodes: A list of nodes on a surface (nnodes, ndim)
         */

        if (nodes.size() < 3){
            std::cerr << "Error: number of nodes must be at least three\n";
            assert(1==0);
        }

        //Compute the centroid of the points on the face
        vectorType faceCentroid;
        vectorTools::computeMean(nodes, faceCentroid);

        std::vector< matrixType > tets;
        tets.reserve(nodes.size());

        unsigned int i=0;
        for (auto it=nodes.begin(); it!=nodes.end()-1; it++, i++){
            tets.push_back({p, faceCentroid, nodes[i], nodes[i+1]});
        }

        tets.push_back({p, faceCentroid, nodes[nodes.size()-1], nodes[0]});

        return tets;
    }

    floatType getTetVolume(const matrixType &tet){
        /*!
         * Compute and return the volume of a tetrahedron
         * 
         * :param const matrixType &tet: The nodes of the tetrahedron
         */

        if (tet.size() != 4){
            std::cerr << "Error: A tetrahedron must be defined by 4 nodes not " << tet.size() << "\n";
            assert(1==0);
        }

        //Define the tetrahedron's sides
        vectorType s1 = tet[1] - tet[0];
        vectorType s2 = tet[2] - tet[0];
        vectorType s3 = tet[3] - tet[0];

        //Compute the volume
        return std::abs(vectorTools::dot(vectorTools::cross(s1, s2), s3))/6;
    }

    int getUnitToTetMap(const matrixType &nodes, matrixType &A, vectorType &d){
        /*!
         * Get the map from the unit tetrahedron to an arbitrary
         * tetrahedron with the given nodes.
         * 
         * :param const matrixType &nodes: The nodes of the tetrahedron 
         *     (nnoddes, ndim)
         * :param matrixType &A: The mapping matrix
         * :param vectorType &d: The shift vector
         */

        if (nodes.size() != 4){
            std::cerr << "Error: A tetrahedron must be defined by 4 nodes not " << nodes.size() << "\n";
            assert(1==0);
        }

        d = nodes[0];

        if (d.size() != 3){
            std::cerr << "Error: A tetrahedron must be defined in 3D not " << d.size() << "\n";
            assert(1==0);
        }

        A.resize(nodes[0].size(), {0, 0, 0});

        for (unsigned int i=1; i<4; i++){
            for (unsigned int j=0; j<3; j++){
                A[j][i-1] = nodes[i][j] - d[j];
            }
        }
        return 0;
    }

    int getTetQuadrature(const unsigned int order, matrixType &points, vectorType &weights){
        /*!
         * Get the local quadrature points and weights of the given order for a tetrahedron
         * Obtained from: https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_mapping_tetrahedrons.pdf
         * 
         * :param int order: The order of integration (0-3 supported)
         */

        if (order==0){
            points = {{0.333333, 0.333333, 0.333333}};
            weights = {1};
        }
        else if (order==1){
            points = {{0.5854101966249685, 0.1381966011250105, 0.1381966011250105},
                      {0.1381966011250105, 0.1381966011250105, 0.1381966011250105},
                      {0.1381966011250105, 0.1381966011250105, 0.5854101966249685},
                      {0.1381966011250105, 0.5854101966249685, 0.1381966011250105}};

            weights = {0.25, 0.25, 0.25, 0.25};
        }
        else if (order==2){
            points = {{0.2500000000000000, 0.2500000000000000, 0.2500000000000000},
                      {0.5000000000000000, 0.1666666666666667, 0.1666666666666667},
                      {0.1666666666666667, 0.1666666666666667, 0.1666666666666667},
                      {0.1666666666666667, 0.1666666666666667, 0.5000000000000000},
                      {0.1666666666666667, 0.5000000000000000, 0.1666666666666667}};
            weights = {-.8, 0.45, 0.45, 0.45, 0.45};
        }
        else if (order==3){
            points = {{0.5684305841968444, 0.1438564719343852, 0.1438564719343852},\
                      {0.1438564719343852, 0.1438564719343852, 0.1438564719343852},\
                      {0.1438564719343852, 0.1438564719343852, 0.5684305841968444},\
                      {0.1438564719343852, 0.5684305841968444, 0.1438564719343852},\
                      {0.0000000000000000, 0.5000000000000000, 0.5000000000000000},\
                      {0.5000000000000000, 0.0000000000000000, 0.5000000000000000},\
                      {0.5000000000000000, 0.5000000000000000, 0.0000000000000000},\
                      {0.5000000000000000, 0.0000000000000000, 0.0000000000000000},\
                      {0.0000000000000000, 0.5000000000000000, 0.0000000000000000},\
                      {0.0000000000000000, 0.0000000000000000, 0.5000000000000000}};
            weights = {0.2177650698804054, 0.2177650698804054, 0.2177650698804054,
                       0.2177650698804054, 0.0214899534130631, 0.0214899534130631,
                       0.0214899534130631, 0.0214899534130631, 0.0214899534130631,
                       0.0214899534130631};
        }
        else{
            std::cerr << "Error: order not supported\n";
            assert(1==0);
        }
        return 0;
    }

    int findPointsOnFace(const vectorType &n, const vectorType &q, const matrixType &points,
                         std::vector< unsigned int > &surfacePoints,
                         double tolr, double tola){
        /*!
         * Determine which points lie on the face
         * 
         * :param const vectorType &n: The face normal
         * :param const vectorType &q: A point on the face
         * :param const matrixType &points: The points to search (npoints, ndim)
         * :param double tolr: The relative tolerance
         * :param double tola: The absolute tolerance
         */

        //Set the tolerance
        double tol = tolr*vectorTools::l2norm(q) + tola;
        
        //Clear out any data in surfacePoints
        surfacePoints.clear();

        //Begin iteration through the points
        unsigned int i=0;
        double d;
        for (auto it = points.begin(); it!=points.end(); it++, i++){
            d = vectorTools::dot(n, *it - q);
            if (abs(d)<=tol){
                surfacePoints.push_back(i);
            }
        }
        return 0;
    }

    int orderPlanarPoints(const matrixType &points, std::vector< unsigned int > &orderedIndices){
        /*!
         * Order a collection of points such that they are ordered in a CCW fashion.
         * 
         * :param const matrixType &points: The points to order (npoints, ndim)
         */

        //Make sure there are points in points
        if (points.size()==0){
            orderedIndices = {};
            return 0;
        }
        else if (points.size()==1){
            orderedIndices = {0};
            return 0;
        }
        else if (points.size()==2){
            orderedIndices = {0, 1};
            return 0;
        }

        //Compute the centroid
        vectorType c = vectorTools::computeMean(points);

        //Compute the unit vector pointing to node 1
        vectorType d = points[0] - c;
        d /= vectorTools::l2norm(d);

        //Compute the unit normal vector
        vectorType n = vectorTools::cross(d, points[1] - c);
        if (vectorTools::fuzzyEquals(vectorTools::l2norm(n), 0.)){
            std::cerr << "Error: the normal vector has zero length. The points are not co-planar.\n";
            assert(1==0);
        }
        n /= vectorTools::l2norm(n);
    
        //Compute the unit orthogonal vector to node 1
        vectorType e = vectorTools::cross(n, d);
        e /= vectorTools::l2norm(e);

        //angle vector
        vectorType angles;
        angles.reserve(points.size());
        angles.push_back(0);

        //Iterate through the remaining points
        vectorType f;
        floatType x, y;
        for (auto p = points.begin() + 1; p!=points.end(); p++){
            f = *p - c;
            x = vectorTools::dot(d, f);
            y = vectorTools::dot(e, f);
            angles.push_back(std::atan2(y, x));
        }

        //Get the indices required to sort the array
        orderedIndices = vectorTools::argsort(angles);
        return 0;
    }

    int getFacePoints(const std::vector< faceType > &faces, const matrixType &points,
        std::vector< std::vector< unsigned int > > &indexFaces){
        /*!
         * Collect the points located on each face ordered CCW
         * 
         * :param const std::vector< faceType > &faces: A vector of faces defined by pairs 
         *     of vectorTypes (normal, point) where normal is the normal of the face and 
         *     and point is a point on the face.
         * :param matrixType points: A collection of points (npoints, ndim)
         * :param std::vector< std::vector< unsigned int > > &indexFaces: A vector of the indices 
         *     of the points located on each face (nfaces, npointsOnFace).
         */

        //Reserve enough memory for the indices
        indexFaces.resize(faces.size());

        std::vector< unsigned int > pointsOnFace;
        std::vector< unsigned int > argSortPoints;
        std::vector< unsigned int > orderedPointsOnFace;
        matrixType subPoints;

        //Iterate through the faces
        unsigned int i=0;
        for (auto face=faces.begin(); face!=faces.end(); face++, i++){
            findPointsOnFace(face->first, face->second, points, pointsOnFace);

            if (pointsOnFace.size() <=3){ //We don't need to sort the points if there are less than four
                indexFaces[i] = pointsOnFace;
            }
            else{
                //Get the points on the face
                vectorTools::getValuesByIndex(points, pointsOnFace, subPoints);
                
                //Sort the points CCW and store them
                orderPlanarPoints(subPoints, argSortPoints);
                vectorTools::getValuesByIndex(pointsOnFace, argSortPoints, orderedPointsOnFace);
                indexFaces[i] = orderedPointsOnFace;
            }
        }
        return 0;
    }

    int volumeToTets(const std::vector< faceType > &faces, const matrixType &points,
        std::vector< matrixType > &tets){
        /*!
         * Deconstruct a volume into a collection of tetrahedra
         * 
         * :param std::vector< faceType > &faces: A vector of faces defined by pairs 
         *     of vectorTypes (normal, point) where normal is the normal of the face and 
         *     and point is a point on the face.
         * :param matrixType points: A collection of points (npoints, ndim)
         * :param std::vector< matrixType > &tets: The tetrahedra which describe the volume
         */

        //Compute the centroid of the points
        vectorType c = vectorTools::computeMean(points);

        //Get the indices of the points associated with each face
        std::vector< std::vector< unsigned int > > facePointIndices;
        getFacePoints(faces, points, facePointIndices);

        //Compute the tetrahedra
        tets.clear();
        matrixType facePoints;
        std::vector< matrixType > faceTets;

        for (auto fpi=facePointIndices.begin(); fpi!=facePointIndices.end(); fpi++){
            if ((*fpi).size() < 3){
                continue;
            }
            vectorTools::getValuesByIndex(points, *fpi, facePoints);
            faceTets = getTets(c, facePoints);

            unsigned int i0 = tets.size();            
            tets.resize(tets.size() + faceTets.size());
            for (auto fT=faceTets.begin(); fT!=faceTets.end(); fT++, i0++){
                tets[i0] = *fT;
            }
        }
        return 0;
    }

    int findMidpoints(const vectorType &p, const matrixType &points, matrixType &midpoints, 
        floatType tolr, floatType tola){
        /*!
         * Find the midpoints between the point p and a collection of points removing any 
         * midpoints which have a distance of zero away from the point p
         * 
         * :param const vectorType &p: The origin point (ndim)
         * :param const matrixType &points: The points to interpolate betwen (npoints, ndim)
         * :param matrixType &midpoints: The midpoints
         * :param floatType tolr: The relative tolerance
         * :param floatType tola: The absolute tolerance
         */

        //Compute the distances
        vectorType distances(points.size(), 0);
        floatType meanDistance = 0;

        unsigned int i=0;
        for (auto point=points.begin(); point!=points.end(); point++, i++){
            distances[i] = vectorTools::l2norm(*point - p);
            meanDistance += distances[i]/points.size();
        }

        //Compute the tolerance
        floatType tol = tolr*meanDistance + tola;

        //Compute the midpoints
        i=0;
        midpoints.clear();
        for (auto point=points.begin(); point!=points.end(); point++, i++){
            if (distances[i]>=tol){
                midpoints.push_back(0.5*((*point) - p) + p);
            }
        }
        return 0;
    }

    int findPointOfIntersection(const std::vector< faceType > &planes, vectorType &point, bool &solveFlag){
        /*!
         * Find the point of intersection of three planes if it exists
         * 
         * :param const std::vector< faceType > &planes: A vector of pairs of the form (normal, point)
         *     where normal is the normal of the surface and point is a point on that surface
         * :param vectorType &point: The point of intersection
         * :param bool &solveFlag: A flag indicating that a unique solution was found
         */

        if (planes.size() != 3){
            std::cerr << "Error: Three planes must be provided\n";
            assert(1==0);
            return 1;
        }

        //Form the matrix equations
        vectorType Avec;
        vectorType b(3);

        unsigned int i=0;
        for (auto plane=planes.begin(); plane!=planes.end(); plane++, i++){
            Avec.insert(Avec.end(), plane->first.begin(), plane->first.end());
            b[i] = vectorTools::dot(plane->first, plane->second);
        }

        unsigned int rank;
        point = vectorTools::solveLinearSystem(Avec, b, 3, 3, rank);
        
        if (rank != 3){
            solveFlag = false;
        }
        else{
            solveFlag = true;
        }
        return 0;
    }

    int findAllPointsOfIntersection(const std::vector< faceType > &planes, matrixType &intersectionPoints,
        floatType tolr, floatType tola){
        /*!
         * Find all of the points of intersection of the set of planes removing duplicates
         * 
         * :param const std::vector< faceType > &planes: A vector of planes of the form (normal, point)
         *     where normal is the normal of the surface and point is a point on that surface.
         * :param matrixType &intersectoinPoints: The points of intersection between the faces.
         * :param floatType tolr: The relative tolerance
         * :param floatType tola: The absolute tolerance
         */

        intersectionPoints.clear();
        vectorType point;
        bool intersects;

        unsigned int i=0;
        for (auto p=planes.begin(); p!=planes.end()-2; p++, i++){

            unsigned int j=i+1;
            for (auto q=planes.begin()+i+1; q!=planes.end()-1; q++, j++){

                if (!vectorTools::isParallel(p->first, q->first)){

                    unsigned int k=j+1;
                    for (auto r=planes.begin()+j+1; r!=planes.end(); r++, k++){

                        findPointOfIntersection({*p, *q, *r}, point, intersects);

                        if (intersects){

                            if (!isDuplicate(point, intersectionPoints)){
                                intersectionPoints.resize(intersectionPoints.size() + 1);
                                intersectionPoints.back() = point;
                            }
                        }
                    }
                }
            }
        }
        return 0;
    }

    bool isDuplicate(const vectorType &point, const matrixType &points){
        /*!
         * Check to see if point is a duplicate of a value in points
         * 
         * :param const vectorType &point: The point in question
         * :param const matrixType &points: The collection of points to check
         */

        bool duplicate=false;

        for (auto p=points.begin(); p!=points.end(); p++){
            duplicate = vectorTools::fuzzyEquals(*p, point);
            if (duplicate){return duplicate;}
        }
        return duplicate;
    }

    int determineInteriorPoints(const vectorType &pInside, const matrixType &points, const std::vector< faceType > &faces,
        std::vector< unsigned int > &interiorPoints, floatType tolr, floatType tola){
        /*!
         * Determine which points are interior points
         * 
         * :param const vectorType &pInside: A point which is known to be located inside of the body
         * :param const matrixType &points: The points to test (npoints, ndim)
         * :param const std::vector< faceType > &faces: A vector of planes of the form (normal, point)
         *     where normal is the normal of the surface and point is a point on that surface
         * :param std::vector< unsigned int > &interiorPoints: The indices of the points which have 
         *     been identified as being inside or on the surface of the domain.
         * :param floatType tolr: The relative tolerance
         * :param floatType tola: The absolute tolerance
         */

        interiorPoints.clear();
        vectorType d, e;
        bool isInside = true;
        double tol;

        unsigned int i=0;
        for (auto p=points.begin(); p!=points.end(); p++, i++){
            d = *p - pInside;
            isInside = true;

            for (auto face=faces.begin(); face!=faces.end(); face++){
                e = face->second - pInside;

                tol = std::fmax(vectorTools::l2norm(d), vectorTools::l2norm(e))*tolr + tola;

                if (tol<vectorTools::dot(face->first, d-e)){
                    isInside = false;
                    break;
                }
            }

            if (isInside){
                interiorPoints.resize(interiorPoints.size() + 1);
                interiorPoints.back() = i;
            }
        }
        return 0;
    }

    int midpointsToFaces(const vectorType &p, const matrixType &midpoints, std::vector< faceType > &faces){
        /*!
         * Convert a collection of midpoints between a set of points to
         * a collection of faces relative to point p
         * 
         * :param const vectorType &p: The point which is located inside of the faces (ndim)
         * :param const matrixType &midpoints: The points in question (npoints, ndim)
         * :param std::vector< faceType > &faces: The faces which bound the domain.
         */

        //Initialize values
        vectorType normal;
        vectorType point;

        faces.clear();

        for (auto mp=midpoints.begin(); mp!=midpoints.end(); mp++){
            normal = (*mp - p)/vectorTools::l2norm(*mp - p);
            faces.push_back(faceType(normal, *mp));
        }
        return 0;
    }

    int getVolumeSubdomainAsTets(const unsigned int index, const matrixType &domainPoints, const std::vector< faceType > &faces,
        std::vector< matrixType > &subdomainTets){
        /*!
         * Get a subdomain of the volume in the form of tetrahedra
         * 
         * :param const unsigned int index: The index of the domain points to use to compute the subdomain.
         * :param const matrixType &domainPoints: The collection of points within each of the subdomains.
         * :param std::vector< faceType > &faces: A vector of planes of the form (normal, point)
         *     where normal is the normal of the surface and point is a point on that surface
         * :param std::vector< matrixType > &subdomainTets: The subdomain broken up into tetrahedra
         */

        //Clear the subdomain
        subdomainTets.clear();

        //Find the midpoints between the domain points
        matrixType midpoints;
        findMidpoints(domainPoints[index], domainPoints, midpoints);

        //Find the faces corresponding to the midpoints and join them with the faces
        std::vector< faceType > domainFaces;
        midpointsToFaces(domainPoints[index], midpoints, domainFaces);
        domainFaces.insert(domainFaces.end(), faces.begin(), faces.end());
        removeDuplicateFaces(domainFaces);

        //Find all of the points of intersection
        matrixType extremePoints;
        findAllPointsOfIntersection(domainFaces, extremePoints);

        //Find the subset of the extreme points that are interior points
        std::vector< unsigned int > interiorPointsIndices;
        determineInteriorPoints(domainPoints[index], extremePoints, domainFaces, interiorPointsIndices);

        //Extract the interior points
        matrixType interiorPoints;
        vectorTools::getValuesByIndex(extremePoints, interiorPointsIndices, interiorPoints);

        //Get the tetrahedra of the volume
        volumeToTets(domainFaces, interiorPoints, subdomainTets);
        return 0;
    }

    int mapLocalTetPointsToGlobal(const matrixType &tet, const matrixType &localPoints, matrixType &globalPoints){
        /*!
         * Map points defined in the tetrahedron's local coordinates to their global coordinates
         * 
         * :param const matrixType &tet: The global nodes of the tetrahedron
         * :param const matrixType &localPoints: The local coordinates of the points to be mapped
         * :param matrixType &globalPoints: The global coordinates of the points
         */

        //Resize globalPoints
        globalPoints.resize(localPoints.size());

        //The mapping terms
        matrixType A;
        vectorType d;

        //Get the map
        getUnitToTetMap(tet, A, d);

        //Apply the map to the points
        unsigned int n=0;
        for (auto point=localPoints.begin(); point!=localPoints.end(); point++, n++){
            //Apply the mapping
            globalPoints[n] = vectorTools::dot(A, *point) + d;
        }
        return 0;
    }

    int writeTetsToFile(const std::string &fileName, const std::vector< matrixType > &tets){
        /*!
         * Write the tets to a file. The format is a list of three comma separated coordinates 
         * followed by a carriage return. Every four points defines a tet.
         * 
         * :param const std::string &fileName: The filename to write to
         * :param const std::vector< matrixType > &tets: The tetrahedra to be written to a file.
         */

        std::ofstream outFile(fileName);
        for (auto tet=tets.begin(); tet!=tets.end(); tet++){
            for (auto point=(*tet).begin(); point!=(*tet).end(); point++){
                outFile << (*point)[0];
                for (auto pi=(*point).begin()+1; pi!=(*point).end(); pi++){
                    outFile << ", " << (*pi);
                }
                outFile << "\n";
            }
        }
        outFile.close();
        return 0;
    }

    int readTetsFromFile(const std::string &fileName, std::vector< matrixType > &tets){
        /*!
         * Read a collection of tets from a file. The format is a list of three comma separated coordinates 
         * followed by a carriage return. Every four points defines a tet.
         * 
         * :param const std::string &fileName: The filename to read from
         * :param const std::vector< matrixType > &tets: The tetrahedra read from the file.
         */

        std::ifstream testTets(fileName);
        std::string line;
        std::stringstream ss;
        std::string subStr;
        std::vector< floatType > point;
        matrixType tet(4);

        unsigned int iter=0;
        if (testTets.is_open()){
            while (std::getline(testTets, line)){
                point.clear();
                if (iter>=4){
                    tets.resize(tets.size() + 1);
                    tets.back() = tet;
                    iter = 0;
                }

                ss = std::stringstream(line);

                while (ss.good()){
                    std::getline(ss, subStr, ',' );
                    point.push_back(std::stod(subStr));
                }
                tet[iter] = point; 

                iter++;
            }     

            testTets.close();
            std::remove(fileName.c_str());
        }
        else{
            std::cerr << "Error: file cannot be opened.\n";
            return 1;
        }
        return 0;
    }

    int removeDuplicateFaces(std::vector< faceType > &faces, const floatType tolr, const floatType tola){
        /*!
         * Remove duplicate faces from a vector of faces
         * 
         * :param std::vector< faceType > : A vector of planes of the form (normal, point)
         *     where normal is the normal of the surface and point is a point on that surface
         * :param const floatType tolr: The relative tolerance
         * :param const floatType tola: The absolute tolerance
         */

        std::vector< unsigned int > uniqueFaces;
        uniqueFaces.reserve(faces.size());
        uniqueFaces.push_back(0);

        double tol, d;
        bool isUnique;

        unsigned int i=1;
        for (auto face=faces.begin()+1; face!=faces.end(); face++, i++){
            isUnique = true;
            //Compare the current face to all previous faces
            for (auto compFace=faces.begin(); compFace!=faces.begin()+i; compFace++){
                //Check if the normals are parallel (I'm not going to allow faces with opposite 
                //normals but which are otherwise identical)
                if (vectorTools::isParallel(face->first, compFace->first)){
                    //Check if the distance between the planes is zero
                    tol = tolr*std::max(vectorTools::l2norm(face->second), vectorTools::l2norm(compFace->second)) + tola;
                    
                    d = vectorTools::dot(face->first, compFace->second - face->second);
                    if (abs(d)<=tol){
                        isUnique = false;
                        break;
                    }
                }
            }
            if (isUnique){
                uniqueFaces.push_back(i);
            }
        }

        std::vector< faceType > unique;
        vectorTools::getValuesByIndex(faces, uniqueFaces, unique);
        faces = unique;
        return 1;
    }
}
