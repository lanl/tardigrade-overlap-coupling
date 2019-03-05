/*!
===============================================================================
|                              element_library.h                              |
===============================================================================
| Header file for the element library. This contains definitions for elements |
| and supporting classes which help construct finite element based solutions. |
===============================================================================
*/

#ifndef ELEMENT_LIBRARY_H
#define ELEMENT_LIBRARY_H

#include<iostream>
#include<vector>
#include<math.h>
#include<assert.h>


namespace elementlib{
    //!===
    //! | Function definitions
    //!===

    //template< typename T1, typename T2 >
    bool fuzzy_compare(const double &a, const double &b, double tolr = 1e-9, double tola = 1e-9);

    //!===
    //! | Class definitions
    //!===

    //!Forward declarations
    class Vector;
    class BaseElement;

    //!Type definitions
    typedef std::vector< std::vector< double > > stdCoordinates;
    typedef std::vector< Vector > VectorCoordinates;

    class BaseElement{
        /*!
        The base for the element class. This class implements the basic 
        functions required to compute the weights and other required 
        terms for the construction of the shape function matrices.

        Inheriting classes are required to define:

        > void initialize(): Initialization which populates local_coordinates, global_coordinates, gauss points, and gauss weights
        > double shape_function(const &int, const &Vector): Function which computes the value of the shape function associated with node n at local position Vector
        > Vector grad_shape_function(const &int, const &Vector): Function which computes the value of the gradient of the shape function associated with node n at local position Vector
        */
    
        public:

            //! > Constructors
            BaseElement(){};
            BaseElement(const stdCoordinates &global_nodes);
            BaseElement(const VectorCoordinates &global_nodes);
            BaseElement(const stdCoordinates &global_nodes, const stdCoordinates &reference_nodes);
            BaseElement(const VectorCoordinates &global_nodes, const VectorCoordinates &reference_nodes);

            //! > Virtual methods

            //! User defined function which populates local_coordinates, global_coordinates,
            //! gauss_points, and gauss_weights
            virtual void initialize() = 0;
            
            //! User defined function which computes the value of the shape function associated with node n at local position Vector
            virtual double shape_function(const int &node, const Vector &Position) = 0;

            //! User defined function which computes the value of the gradient shape function w.r.t. the local coordinates associated with node n at local position Vector.
            virtual Vector grad_shape_function(const int &node, const Vector &Position) = 0;

            //! > Methods

            void interpolate(const std::vector< Vector > &nodal_values, const Vector &Position, Vector &result);

            void local_gradient(const std::vector< Vector > &nodal_values, const Vector &Position, std::vector< Vector > &result);

            void compute_dxdxi(const Vector &Position, std::vector< Vector > &result);

            Vector get_local_coordinates(const int &n) const;

            void print() const;

        protected:

            //!The local coordinates of the nodes
            VectorCoordinates local_coordinates;
            
            //!The global coordinates of the nodes
            VectorCoordinates global_coordinates;

            //!The reference coordinates of the nodes
            VectorCoordinates reference_coordinates;

            //!The local coordinates of the gauss points
            VectorCoordinates gauss_points;

            //!The weights of the gauss points
            std::vector< double > gauss_weights;

    };

    class Hex8: public BaseElement{
        /*!
        A fully integrated linear hexahedron element.
        */

        public:

            Hex8(const stdCoordinates &global_coordinates): BaseElement(global_coordinates){initialize();}
            Hex8(const VectorCoordinates &global_coordinates): BaseElement(global_coordinates){initialize();}
            Hex8(const stdCoordinates &global_coordinates, const stdCoordinates &reference_coordinates): BaseElement(global_coordinates, reference_coordinates){initialize();}
            Hex8(const VectorCoordinates &global_coordinates, const VectorCoordinates &reference_coordinates): BaseElement(global_coordinates, reference_coordinates){initialize();}

            void initialize() override;

            double shape_function(const int &node, const Vector &Position) override;

            Vector grad_shape_function(const int &node, const Vector &Position) override;

    };

    class Vector{
        /*!

        Utility class providing vector addition, subtraction, etc.

        Converts all numbers to double.

        */

        public:

            Vector(){}

            template< typename T >
            Vector(const std::vector< T > &vec){
                _value = vec;
            }

            Vector(const Vector &v){
                _value = v._value;
            }

            void print() const{
                /*!
                Print the vector values to the terminal.
                */

                for (unsigned int i=0; i<_value.size(); i++){
                    std::cout << _value[i] << " ";
                }
                std::cout << "\n";
            }

            bool operator ==(const Vector &b) const{
                /*!
                Equality operator.
                */

                if (_value.size() != b._value.size()){return false;}

                for (unsigned int i=0; i<_value.size(); i++){
                    if (!fuzzy_compare(_value[i], b._value[i])){return false;}
                }

                return true;
            }

            double operator ()(const unsigned int i) const{
                /*!
                Component access operator.
                */

                return _value[i];
            }

            Vector& operator +=(const double &a){
                /*!
                Compound additive assignment
                */
                for (unsigned int i=0; i<this->_value.size(); i++){
                    this->_value[i] += a;
                }
                return *this;
            }

            Vector& operator += (const Vector &b){
                /*!
                Compound additive assignment
                */
                if (this->_value.size() != b._value.size()){
                    std::cout << "Error: Vectors of different sizes cannot be added\n";
                    assert(1==0);
                }

                for (unsigned int i=0; i<this->_value.size(); i++){
                    this->_value[i] += b._value[i];
                }
                return *this;
            }

            Vector operator -() const{
                /*!
                Negation operator
                */

                std::vector< double > vec(_value.size(), 0);
                for (unsigned int i=0; i<_value.size(); i++){
                    vec[i] = -_value[i];
                }
                Vector result(vec);
                return result;
            }

            Vector& operator -=(const Vector &b){
                /*!
                Compound subtractive assignment
                */

                *this += -b;
                return *this;

            }

            Vector& operator -=(const double &b){
                /*!
                Compound subtractive assignment
                */
                *this += -b;
                return *this;

            }

            Vector& operator *=(const double &a){
                /*!
                Compund multiplicative assignment
                */
                for (unsigned int i=0; i<this->_value.size(); i++){
                    this->_value[i] *= a;
                }
                return *this;
            }

            Vector& operator *=(const Vector &a){
                /*!
                Multiplication operator (elementwise product)
                */

                if (this->_value.size() !=a._value.size()){
                    std::cout << "Error: Incompatible shapes\n";
                    assert(1==0);
                }

                for (unsigned int i=0; i<this->_value.size(); i++){
                    this->_value[i] *= a._value[i];
                }
                return *this;
            }

            Vector& operator /=(const double &a){
                /*!
                Compound division assignment.
                */

                *this *=(1./a);
                return *this;
            }

            Vector& operator /=(const int &a){
                /*!
                Compound division assignment.
                */

                *this *= (1./a);
                return *this;
            }

            double sum() const{
                /*!
                Sum the elements of the vector together
                */

                double result = 0;
                for (unsigned int i=0; i<_value.size(); i++){
                    result += _value[i];
                }
                return result;
            }

            double product() const{
                /*!
                Multiply the elements of the vector together
                */

                double result = 1;
                for (unsigned int i=0; i<_value.size(); i++){
                    result *= _value[i];
                }
                return result;
            }

            std::vector< Vector > dyadic_product(const Vector &b) const{
                /*!
                Compute the dyadic product of this vector with the provided vector. 
                */
            
                std::vector< Vector > result(_value.size());
                for (unsigned int i=0; i<this->_value.size(); i++){
                    result[i] = b;
                    result[i] *= this->_value[i];
                }
                return result;
            }

        protected:
            std::vector< double > _value;
        
    };

    //!Vector overloaded operators
    inline Vector operator +(Vector a, const Vector &b){
        /*!
        Addition operator.
        */

        a += b;
        return a;
    }

    inline Vector operator +(const double &a, Vector b){
        /*!
        Addition operator (elementwise addition)
        */

        b += a;
        return b;
    }

    inline Vector operator +(Vector b, const double &a){
        /*!
        Addition operator (elementwise addition)
        */

        b += a;
        return b;
    }

    inline Vector operator -(Vector a, const Vector &b){
        /*!
        Subtraction operator
        */

        a -= b;
        return a;
    }

    inline Vector operator -(Vector a, const double &b){
        /*!
        Subtraction operator
        */

        a -= b;
        return a;
    }

    inline Vector operator -(const double &a, const Vector &b){
        /*!
        Subtraction operator
        */

        Vector c = -b;
        c += a;
        return c;
    }

    inline Vector operator *(Vector a, const double &b){
        /*!
        Multiplication operator
        */
        a *= b;
        return a;
    }

    inline Vector operator *(const double &b, Vector a){
        /*!
        Multiplication operator
        */
        a *= b;
        return a;
    }

    inline Vector operator *(Vector a, const Vector &b){
        /*!
        Multiplication operator (elementwise product)
        */
        a *= b;
        return a;
    }

    inline Vector operator /(Vector a, const double b){
        /*!
        Elementwise division by scalar
        */

        a /= b;
        return a;
    }
}

#endif 
