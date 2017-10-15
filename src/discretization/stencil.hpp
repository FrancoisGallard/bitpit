/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __BTPIT_STENCIL_HPP__
#define __BTPIT_STENCIL_HPP__

#include <array>
#include <unordered_map>

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"

// Stream operators for the base class
namespace bitpit {

template<typename weight_t>
class BaseStencil;

}

template<typename weight_t>
bitpit::OBinaryStream & operator<<(bitpit::OBinaryStream &buffer, const bitpit::BaseStencil<weight_t> &stencil);

template<typename weight_t>
bitpit::IBinaryStream & operator>>(bitpit::IBinaryStream &buffer, bitpit::BaseStencil<weight_t> &stencil);

// Declaration of the base class
namespace bitpit {

template <typename weight_t>
class BaseStencil {

template<typename U>
friend bitpit::OBinaryStream & (::operator<<) (bitpit::OBinaryStream &buffer, const BaseStencil<U> &stencil);
template<typename U>
friend bitpit::IBinaryStream & (::operator>>) (bitpit::IBinaryStream &buffer, BaseStencil<U> &stencil);

public:
    long NULL_ID = - std::numeric_limits<long>::max();

    typedef weight_t weight_type;

    struct item_type {
        long id;
        weight_type weight;
    };

    BaseStencil(const weight_t &zero = weight_t());
    BaseStencil(int nBuckets, const weight_t &zero = weight_t());
    BaseStencil(int nBuckets, int nBucketItems, const weight_t &zero = weight_t());
    BaseStencil(const std::vector<int> &sizes, const weight_t &zero = weight_t());

    void initialize(const weight_t &zero = weight_t());
    void initialize(int nBuckets, const weight_t &zero = weight_t());
    void initialize(int nBuckets, int nBucketItems, const weight_t &zero = weight_t());
    void initialize(const std::vector<int> &bucketSizes, const weight_t &zero = weight_t());

    void clear();

    std::size_t size() const;
    std::size_t size(int bucket) const;

    int getBucketCount() const;

    long & getPattern(std::size_t pos);
    const long & getPattern(std::size_t pos) const;
    long & getPattern(int bucket, std::size_t pos);
    const long & getPattern(int bucket, std::size_t pos) const;
    long * patternData();
    const long * patternData() const;
    const FlatVector2D<long> & getPattern() const;
    void setPattern(std::size_t pos, long id);
    void setPattern(int bucket, std::size_t pos, long id);

    weight_t & getWeight(std::size_t pos);
    const weight_t & getWeight(std::size_t pos) const;
    weight_t & getWeight(int bucket, std::size_t pos);
    const weight_t & getWeight(int bucket, std::size_t pos) const;
    weight_t * weightData();
    const weight_t * weightData() const;
    const FlatVector2D<weight_t> & getWeights() const;
    void setWeight(std::size_t pos, const weight_t &weight);
    void setWeight(int bucket, std::size_t pos, const weight_t &weight);

    item_type getItem(std::size_t pos) const;
    item_type getItem(int bucket, std::size_t pos) const;
    void setItem(std::size_t pos, long id, const weight_t &weight);
    void setItem(int bucket, std::size_t pos, long id, const weight_t &weight);
    void appendItem(long id, const weight_t &weight);
    void appendItem(int bucket, long id, const weight_t &weight);

    weight_t & getConstant();
    const weight_t & getConstant() const;
    void setConstant(const weight_t &value);
    void sumConstant(const weight_t &value);

    void flatten();
    void optimize(double tolerance = 1.e-12);
    void renumber(const std::unordered_map<long, long> &map);
    void addComplementToZero(const long id);

    void display(std::ostream &out, double factor = 1.) const;

    size_t getBinarySize() const;

    BaseStencil<weight_t> & operator*=(double factor);
    BaseStencil<weight_t> & operator/=(double factor);
    BaseStencil<weight_t> & operator+=(const BaseStencil<weight_t> &other);
    BaseStencil<weight_t> & operator-=(const BaseStencil<weight_t> &other);

private:
    weight_t m_zero;
    FlatVector2D<long> m_pattern;
    FlatVector2D<weight_t> m_weights;
    weight_t m_constant;

    weight_t * find(int bucket, long id);
    const weight_t * find(int bucket, long id) const;

    template<typename U = weight_t, typename std::enable_if<std::is_fundamental<U>::value>::type* = nullptr>
    bool isWeightNeglibile(int bucket, size_t k, double tolerance = 1.e-12);

    template<typename U = weight_t, typename std::enable_if<!std::is_fundamental<U>::value>::type* = nullptr>
    bool isWeightNeglibile(int bucket, size_t k, double tolerance = 1.e-12);

};

}

// Operators for the base class
template <typename weight_t>
bitpit::BaseStencil<weight_t> operator*(const bitpit::BaseStencil<weight_t> &stencil, double factor);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator*(double factor, const bitpit::BaseStencil<weight_t> &stencil);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator/(const bitpit::BaseStencil<weight_t> &stencil, double factor);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator+(const bitpit::BaseStencil<weight_t> &, const bitpit::BaseStencil<weight_t> &);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator-(const bitpit::BaseStencil<weight_t> &, const bitpit::BaseStencil<weight_t> &);

// //Stream operators for the specializations
namespace bitpit {

class StencilScalar;
class StencilVector;

}

bitpit::OBinaryStream & operator<<(bitpit::OBinaryStream &buffer, const bitpit::StencilScalar &stencil);
bitpit::IBinaryStream & operator>>(bitpit::IBinaryStream &buffer, bitpit::StencilScalar &stencil);

bitpit::OBinaryStream & operator<<(bitpit::OBinaryStream &buffer, const bitpit::StencilVector &stencil);
bitpit::IBinaryStream & operator>>(bitpit::IBinaryStream &buffer, bitpit::StencilVector &stencil);

// Template implementation
#include "stencil.tpp"

// Declaration of the specializations
namespace bitpit {

class StencilScalar : public BaseStencil<double> {

friend bitpit::OBinaryStream & (::operator<<) (bitpit::OBinaryStream &buffer, const StencilScalar &stencil);
friend bitpit::IBinaryStream & (::operator>>) (bitpit::IBinaryStream &buffer, StencilScalar &stencil);

public:
    typedef BaseStencil::weight_type weight_type;

    StencilScalar();
    StencilScalar(const BaseStencil &);
    StencilScalar(BaseStencil &&);

};

class StencilVector : public BaseStencil<std::array<double, 3>> {

friend bitpit::OBinaryStream & (::operator<<) (bitpit::OBinaryStream &buffer, const StencilVector &stencil);
friend bitpit::IBinaryStream & (::operator>>) (bitpit::IBinaryStream &buffer, StencilVector &stencil);

public:
    typedef BaseStencil::weight_type weight_type;

    StencilVector();
    StencilVector(const BaseStencil &other);
    StencilVector(BaseStencil &&other);

};

}

// Operators for the specializations
bitpit::StencilVector operator*(const bitpit::StencilScalar &, const std::array<double,3> &);

// Functions for the specializations
bitpit::StencilScalar dotProduct(const bitpit::StencilVector &stencil_A, const bitpit::StencilVector::weight_type &vector);

#endif
