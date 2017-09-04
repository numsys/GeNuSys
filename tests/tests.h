/*
GeNuSys - computations with generalized number systems
Copyright (C) 2015-2017  Bence Németh
Copyright (C) 2017  Tamás Krutki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "test_suite.h"
#include "test_utils.h"

#include <GeNuSys/element_traits.h>

#include <GeNuSys/vector.h>
#include <GeNuSys/sparse_vector.h>
#include <GeNuSys/matrix.h>
#include <GeNuSys/sparse_matrix.h>

#include <GeNuSys/p_norm.h>
#include <GeNuSys/frobenius_norm.h>
#include <GeNuSys/operator_norm.h>

class VectorTest : public GeNuSys::Tests::TestSuite
{

    public:

        VectorTest(): TestSuite("Vector") {}

        void run()
        {
            GeNuSys::LinAlg::Vector<int> vector(3);
            assertEqual<unsigned int>(3, vector.getLength(), "Vector length");
            assertEqual(0, vector[0], "Vector component initialization");

            vector.set(2, 1);
            assertEqual(1, vector[2], "Vector component set and access");

            GeNuSys::LinAlg::Vector<int> vectorCopy = vector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(vectorCopy, vector), "Vector copy matches original");

            GeNuSys::LinAlg::Vector<int> largerVector(5);
            largerVector.set(0, 1);
            assertEqual<unsigned int>(5, largerVector.getLength(), "Longer vector length");

            vector = largerVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(vector, largerVector), "Vectors match after assign (operator =)");

            GeNuSys::LinAlg::Vector<double> doubleVector = vector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(doubleVector, vector), "Vectors match after convert and copy");

            doubleVector = largerVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(doubleVector, largerVector), "Vectors match after convert and assign (operator =)");

            GeNuSys::LinAlg::SparseVector<int> sparseVector(4);
            sparseVector.set(0, 1);
            sparseVector.set(3, 2);

            GeNuSys::LinAlg::SparseVector<int> largerSparseVector(6);
            largerSparseVector.set(0, 1);
            largerSparseVector.set(3, 2);
            largerSparseVector.set(5, 3);

            GeNuSys::LinAlg::Vector<int> vectorB = sparseVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(vectorB, sparseVector), "Dense copy of sparse vector matches original");

            vectorB = largerSparseVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(vectorB, largerSparseVector), "Dense vector matches original sparse vector after assign sparse (operator =)");

            GeNuSys::LinAlg::Vector<double> doubleVectorB = sparseVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(doubleVectorB, sparseVector), "Converted dense copy of sparse vector matches original");

            doubleVectorB = largerSparseVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(doubleVectorB, largerSparseVector), "Converted dense vector matches original sparse vector after assign (operator =)");
        }

};

class SparseVectorTest : public GeNuSys::Tests::TestSuite
{

    public:

        SparseVectorTest(): TestSuite("SparseVector") {}

        void run()
        {
            GeNuSys::LinAlg::SparseVector<int> sparseVector(3);
            assertEqual<unsigned int>(3, sparseVector.getLength(), "SparseVector length");
            assertEqual(0, sparseVector[0], "SparseVector component initialization");

            sparseVector.set(2, 1);
            assertEqual(1, sparseVector[2], "SparseVector component set and access");

            GeNuSys::LinAlg::SparseVector<int> sparseVectorCopy = sparseVector;
            assertTrue(GeNuSys::Tests::TestUtils::equals(sparseVectorCopy, sparseVector), "Copy of SparseVector matches original");

            GeNuSys::LinAlg::SparseVector<int> largerSparseVector(5);
            largerSparseVector.set(0, 1);

            assertEqual<unsigned int>(5, largerSparseVector.getLength(), "Longer vector length");

            sparseVector = largerSparseVector;

            assertTrue(GeNuSys::Tests::TestUtils::equals(sparseVector, largerSparseVector), "SparseVector matches original after assign (operator =)");

            GeNuSys::LinAlg::SparseVector<double> doubleSparseVector = sparseVector;

            assertEqual(sparseVector.getLength(), doubleSparseVector.getLength(), "Convereted vector length");

            bool convertedEqual = true;
            for (unsigned int i = 0; i < doubleSparseVector.getLength(); ++i)
            {
                convertedEqual = (doubleSparseVector[i] == (double) sparseVector[i]);
            }
            assertTrue(convertedEqual, "Components of converted match");

            doubleSparseVector = largerSparseVector;

            assertEqual<unsigned int>(5, doubleSparseVector.getLength(), "Correct length after convert and assign (operator =)");

            bool assignConvertEqual = true;
            for (unsigned int i = 0; i < doubleSparseVector.getLength(); ++i)
            {
                assignConvertEqual = (doubleSparseVector[i] == (double) largerSparseVector[i]);
            }
            assertTrue(assignConvertEqual, "Components match after convert and assign (operator =)");

            GeNuSys::LinAlg::Vector<int> vector(4);
            vector.set(0, 1);
            vector.set(3, 2);

            GeNuSys::LinAlg::Vector<int> largerVector(6);
            largerVector.set(0, 1);
            largerVector.set(3, 2);
            largerVector.set(5, 3);

            GeNuSys::LinAlg::SparseVector<int> sparseVectorB = vector;

            assertEqual(vector.getLength(), sparseVectorB.getLength(), "Correct length after copy dense vector");

            bool sparseEqual = true;
            for (unsigned int i = 0; i < sparseVectorB.getLength(); ++i)
            {
                sparseEqual = (sparseVectorB[i] == vector[i]);
            }
            assertTrue(sparseEqual, "Components of dense vector match");

            sparseVectorB = largerVector;

            assertEqual<unsigned int>(largerVector.getLength(), sparseVectorB.getLength(), "Correct length after assign dense vector");

            bool assignSparseEqual = true;
            for (unsigned int i = 0; i < sparseVectorB.getLength(); ++i)
            {
                assignSparseEqual = (sparseVectorB[i] == largerVector[i]);
            }
            assertTrue(assignSparseEqual, "Components match after assign sparse (operator =)");

            GeNuSys::LinAlg::SparseVector<double> doubleSparseVectorB = vector;

            assertEqual(vector.getLength(), doubleSparseVectorB.getLength(), "Correct length after copy and convert dense vector");

            bool sparseConvertEqual = true;
            for (unsigned int i = 0; i < doubleSparseVectorB.getLength(); ++i)
            {
                sparseConvertEqual = (doubleSparseVectorB[i] == (double) vector[i]);
            }
            assertTrue(sparseConvertEqual, "Components of converted dense vector match");

            doubleSparseVectorB = largerVector;

            assertEqual<unsigned int>(largerVector.getLength(), doubleSparseVectorB.getLength(), "Correct length after convert and assign dense");

            bool assignSparseConvertEqual = true;
            for (unsigned int i = 0; i < doubleSparseVectorB.getLength(); ++i)
            {
                assignSparseConvertEqual = (doubleSparseVectorB[i] == (double) largerVector[i]);
            }
            assertTrue(assignSparseConvertEqual, "Components match after convert and assign dense (operator =)");
        }

};

class VectorCmpTest : public GeNuSys::Tests::TestSuite
{

    public:

        VectorCmpTest(): TestSuite("Vector compare") {}

        void run()
        {
            GeNuSys::LinAlg::Vector<int> vector(3);
            vector.set(1, 1);

            GeNuSys::LinAlg::Vector<int> vectorCopy = vector;

            assertTrue(vector == vectorCopy, "Compare with self (operator ==)");
            assertFalse(vector < vectorCopy, "Compare with self (operator <)");
            assertTrue(vector <= vectorCopy, "Compare with self (operator <=)");
            assertFalse(vector > vectorCopy, "Compare with self (operator >)");
            assertTrue(vector >= vectorCopy, "Compare with self (operator >=)");

            for (unsigned int i = 0; i < vectorCopy.getLength(); ++i)
            {
                vectorCopy.set(i, vectorCopy[i] + 1);
            }

            assertFalse(vector == vectorCopy, "Compare with larger (operator ==)");
            assertTrue(vector < vectorCopy, "Compare with larger (operator <)");
            assertTrue(vector <= vectorCopy, "Compare with larger (operator <=)");
            assertFalse(vector > vectorCopy, "Compare with larger (operator >)");
            assertFalse(vector >= vectorCopy, "Compare with larger (operator >=)");

            for (unsigned int i = 0; i < vectorCopy.getLength(); ++i)
            {
                vectorCopy.set(i, vectorCopy[i] - 2);
            }

            assertFalse(vector == vectorCopy, "Compare with smaller (operator ==)");
            assertFalse(vector < vectorCopy, "Compare with smaller (operator <)");
            assertFalse(vector <= vectorCopy, "Compare with smaller (operator <=)");
            assertTrue(vector > vectorCopy, "Compare with smaller (operator >)");
            assertTrue(vector >= vectorCopy, "Compare with smaller (operator >=)");

            GeNuSys::LinAlg::SparseVector<int> vectorSparseCopy = vector;

            assertTrue(vector == vectorSparseCopy, "Compare with sparse copy (operator ==)");
            assertFalse(vector < vectorSparseCopy, "Compare with sparse copy (operator <)");
            assertTrue(vector <= vectorSparseCopy, "Compare with sparse copy (operator <=)");
            assertFalse(vector > vectorSparseCopy, "Compare with sparse copy (operator >)");
            assertTrue(vector >= vectorSparseCopy, "Compare with sparse copy (operator >=)");

            for (unsigned int i = 0; i < vector.getLength(); ++i)
            {
                vector.set(i, vector[i] - 1);
            }

            assertFalse(vector == vectorSparseCopy, "Compare with larger sparse (operator ==)");
            assertTrue(vector < vectorSparseCopy, "Compare with larger sparse (operator <)");
            assertTrue(vector <= vectorSparseCopy, "Compare with larger sparse (operator <=)");
            assertFalse(vector > vectorSparseCopy, "Compare with larger sparse (operator >)");
            assertFalse(vector >= vectorSparseCopy, "Compare with larger sparse (operator >=)");

            for (unsigned int i = 0; i < vector.getLength(); ++i)
            {
                vector.set(i, vector[i] + 2);
            }

            assertFalse(vector == vectorSparseCopy, "Compare with smaller sparse (operator ==)");
            assertFalse(vector < vectorSparseCopy, "Compare with smaller sparse (operator <)");
            assertFalse(vector <= vectorSparseCopy, "Compare with smaller sparse (operator <=)");
            assertTrue(vector > vectorSparseCopy, "Compare with smaller sparse (operator >)");
            assertTrue(vector >= vectorSparseCopy, "Compare with smaller sparse (operator >=)");
        }

};

class SparseVectorCmpTest : public GeNuSys::Tests::TestSuite
{

    public:

        SparseVectorCmpTest(): TestSuite("SparseVector compare") {}

        void run()
        {
            GeNuSys::LinAlg::SparseVector<int> vector(3);
            vector.set(1, 1);

            GeNuSys::LinAlg::SparseVector<int> vectorCopy = vector;

            assertTrue(vector == vectorCopy, "Compare with self (operator ==)");
            assertFalse(vector < vectorCopy, "Compare with self (operator <)");
            assertTrue(vector <= vectorCopy, "Compare with self (operator <=)");
            assertFalse(vector > vectorCopy, "Compare with self (operator >)");
            assertTrue(vector >= vectorCopy, "Compare with self (operator >=)");

            for (unsigned int i = 0; i < vectorCopy.getLength(); ++i)
            {
                vectorCopy.set(i, vectorCopy[i] + 1);
            }

            assertFalse(vector == vectorCopy, "Compare with larger (operator ==)");
            assertTrue(vector < vectorCopy, "Compare with larger (operator <)");
            assertTrue(vector <= vectorCopy, "Compare with larger (operator <=)");
            assertFalse(vector > vectorCopy, "Compare with larger (operator >)");
            assertFalse(vector >= vectorCopy, "Compare with larger (operator >=)");

            for (unsigned int i = 0; i < vectorCopy.getLength(); ++i)
            {
                vectorCopy.set(i, vectorCopy[i] - 2);
            }

            assertFalse(vector == vectorCopy, "Compare with smaller (operator ==)");
            assertFalse(vector < vectorCopy, "Compare with smaller (operator <)");
            assertFalse(vector <= vectorCopy, "Compare with smaller (operator <=)");
            assertTrue(vector > vectorCopy, "Compare with smaller (operator >)");
            assertTrue(vector >= vectorCopy, "Compare with smaller (operator >=)");

            GeNuSys::LinAlg::Vector<int> vectorDenseCopy = vector;

            assertTrue(vector == vectorDenseCopy, "Compare with dense copy (operator ==)");
            assertFalse(vector < vectorDenseCopy, "Compare with dense copy (operator <)");
            assertTrue(vector <= vectorDenseCopy, "Compare with dense copy (operator <=)");
            assertFalse(vector > vectorDenseCopy, "Compare with dense copy (operator >)");
            assertTrue(vector >= vectorDenseCopy, "Compare with dense copy (operator >=)");

            for (unsigned int i = 0; i < vectorDenseCopy.getLength(); ++i)
            {
                vectorDenseCopy.set(i, vectorDenseCopy[i] + 1);
            }

            assertFalse(vector == vectorDenseCopy, "Compare with larger dense (operator ==)");
            assertTrue(vector < vectorDenseCopy, "Compare with larger dense (operator <)");
            assertTrue(vector <= vectorDenseCopy, "Compare with larger dense (operator <=)");
            assertFalse(vector > vectorDenseCopy, "Compare with larger dense (operator >)");
            assertFalse(vector >= vectorDenseCopy, "Compare with larger dense (operator >=)");

            for (unsigned int i = 0; i < vectorDenseCopy.getLength(); ++i)
            {
                vectorDenseCopy.set(i, vectorDenseCopy[i] - 2);
            }

            assertFalse(vector == vectorDenseCopy, "Compare with smaller dense (operator ==)");
            assertFalse(vector < vectorDenseCopy, "Compare with smaller dense (operator <)");
            assertFalse(vector <= vectorDenseCopy, "Compare with smaller dense (operator <=)");
            assertTrue(vector > vectorDenseCopy, "Compare with smaller dense (operator >)");
            assertTrue(vector >= vectorDenseCopy, "Compare with smaller dense (operator >=)");
        }

};

class VectorOperatorTest : public GeNuSys::Tests::TestSuite
{

    public:

        VectorOperatorTest(): TestSuite("Vector Operators") {}

        void run()
        {
            GeNuSys::LinAlg::Vector<int> vector(5);
            for (unsigned int i = 0; i < vector.getLength(); ++i)
            {
                vector.set(i, i + 1);
            }

            GeNuSys::LinAlg::Vector<double> divVector = vector / 2;

            assertEqual(divVector.getLength(), vector.getLength(), "Division doesn't change vector length");

            bool divEquals = true;
            for (unsigned int i = 0; i < divVector.getLength(); ++i)
            {
                divEquals = (divVector[i] == ((double)vector[i] / 2));
            }
            assertTrue(divEquals, "Divide with constant (operator /)");

            GeNuSys::LinAlg::Vector<int> mulVector = vector * 2;

            assertEqual(mulVector.getLength(), vector.getLength(), "Multiplication doesn't change vector length");

            bool mulEquals = true;
            for (unsigned int i = 0; i < mulVector.getLength(); ++i)
            {
                mulEquals = (mulVector[i] == (vector[i] * 2));
            }
            assertTrue(mulEquals, "Multiplication with constant (operator *)");

            GeNuSys::LinAlg::Vector<int> vectorB(5);
            for (unsigned int i = 0; i < vectorB.getLength(); ++i)
            {
                vectorB.set(i, i + 2);
            }

            GeNuSys::LinAlg::Vector<int> addVector = vector + vectorB;

            assertEqual(addVector.getLength(), vector.getLength(), "Add doesn't change vector length");

            bool addEquals = true;
            for (unsigned int i = 0; i < addVector.getLength(); ++i)
            {
                addEquals = (addVector[i] == vector[i] + vectorB[i]);
            }
            assertTrue(addEquals, "Add (operator +)");

            GeNuSys::LinAlg::Vector<int> subVector = vector - vectorB;

            assertEqual(subVector.getLength(), vector.getLength(), "Subtract doesn't change vector length");

            bool subEquals = true;
            for (unsigned int i = 0; i < subVector.getLength(); ++i)
            {
                subEquals = (subVector[i] == vector[i] - vectorB[i]);
            }
            assertTrue(subEquals, "Subtract (operator -)");

            assertEqual(70, vector * vectorB, "Scalar product (operator *)");

            vectorB.set(0, 3);
            assertEqual(71, vector * vectorB, "Scalar product, modified (operator *)");

            GeNuSys::LinAlg::SparseVector<int> sparseVectorB(5);
            for (unsigned int i = 0; i < sparseVectorB.getLength(); ++i)
            {
                if (i % 2 == 0)
                {
                    sparseVectorB.set(i, i + 2);
                }
            }

            GeNuSys::LinAlg::Vector<int> addSparseVector = vector + sparseVectorB;

            assertEqual(addSparseVector.getLength(), vector.getLength(), "Add sparse doesn't change vector length");

            bool addSparseEquals = true;
            for (unsigned int i = 0; i < addSparseVector.getLength(); ++i)
            {
                addSparseEquals = (addSparseVector[i] == vector[i] + sparseVectorB[i]);
            }
            assertTrue(addSparseEquals, "Add sparse (operator +)");

            GeNuSys::LinAlg::Vector<int> subSparseVector = vector - sparseVectorB;

            assertEqual(subSparseVector.getLength(), vector.getLength(), "Subtract sparse doesn't change vector length");

            bool subSparseEquals = true;
            for (unsigned int i = 0; i < subSparseVector.getLength(); ++i)
            {
                subSparseEquals = (subSparseVector[i] == vector[i] - sparseVectorB[i]);
            }
            assertTrue(subSparseEquals, "Subtract sparse (operator -)");

            assertEqual(44, vector * sparseVectorB, "Scalar product (operator *)");

            sparseVectorB.set(0, 3);
            assertEqual(45, vector * sparseVectorB, "Scalar product, modified (operator *)");
        }

};

class SparseVectorOperatorTest : public GeNuSys::Tests::TestSuite
{

    public:

        SparseVectorOperatorTest(): TestSuite("SparseVector Operators") {}

        void run()
        {
            GeNuSys::LinAlg::SparseVector<int> vector(5);
            for (unsigned int i = 0; i < vector.getLength(); ++i)
            {
                if (i % 2 == 0)
                {
                    vector.set(i, i + 1);
                }
            }

            GeNuSys::LinAlg::SparseVector<double> divVector = vector / 2;

            assertEqual(divVector.getLength(), vector.getLength(), "Division doesn't change vector length");

            bool divEquals = true;
            for (unsigned int i = 0; i < divVector.getLength(); ++i)
            {
                divEquals = (divVector[i] == ((double)vector[i] / 2));
            }
            assertTrue(divEquals, "Divide with constant (operator /)");

            GeNuSys::LinAlg::SparseVector<int> mulVector = vector * 2;

            assertEqual(mulVector.getLength(), vector.getLength(), "Multiplication doesn't change vector length");

            bool mulEquals = true;
            for (unsigned int i = 0; i < mulVector.getLength(); ++i)
            {
                mulEquals = (mulVector[i] == (vector[i] * 2));
            }
            assertTrue(mulEquals, "Multiplication with constant (operator *)");

            GeNuSys::LinAlg::Vector<int> vectorB(5);
            for (unsigned int i = 0; i < vectorB.getLength(); ++i)
            {
                vectorB.set(i, i + 2);
            }

            GeNuSys::LinAlg::Vector<int> addVector = vector + vectorB;

            assertEqual(addVector.getLength(), vector.getLength(), "Add doesn't change vector length");

            bool addEquals = true;
            for (unsigned int i = 0; i < addVector.getLength(); ++i)
            {
                addEquals = (addVector[i] == vector[i] + vectorB[i]);
            }
            assertTrue(addEquals, "Add (operator +)");

            GeNuSys::LinAlg::Vector<int> subVector = vector - vectorB;

            assertEqual(subVector.getLength(), vector.getLength(), "Subtract doesn't change vector length");

            bool subEquals = true;
            for (unsigned int i = 0; i < subVector.getLength(); ++i)
            {
                subEquals = (subVector[i] == vector[i] - vectorB[i]);
            }
            assertTrue(subEquals, "Subtract (operator -)");

            assertEqual(44, vector * vectorB, "Scalar product (operator *)");

            vectorB.set(0, 3);
            assertEqual(45, vector * vectorB, "Scalar product, modified (operator *)");

            GeNuSys::LinAlg::SparseVector<int> sparseVectorB(5);
            for (unsigned int i = 0; i < sparseVectorB.getLength(); ++i)
            {
                if (i % 2 == 0)
                {
                    sparseVectorB.set(i, i + 2);
                }
            }

            GeNuSys::LinAlg::SparseVector<int> addSparseVector = vector + sparseVectorB;

            assertEqual(addSparseVector.getLength(), vector.getLength(), "Add sparse doesn't change vector length");

            bool addSparseEquals = true;
            for (unsigned int i = 0; i < addSparseVector.getLength(); ++i)
            {
                addSparseEquals = (addSparseVector[i] == vector[i] + sparseVectorB[i]);
            }
            assertTrue(addSparseEquals, "Add sparse (operator +)");

            GeNuSys::LinAlg::SparseVector<int> subSparseVector = vector - sparseVectorB;

            assertEqual(subSparseVector.getLength(), vector.getLength(), "Subtract sparse doesn't change vector length");

            bool subSparseEquals = true;
            for (unsigned int i = 0; i < subSparseVector.getLength(); ++i)
            {
                subSparseEquals = (subSparseVector[i] == vector[i] - sparseVectorB[i]);
            }
            assertTrue(subSparseEquals, "Subtract sparse (operator -)");

            assertEqual(44, vector * sparseVectorB, "Scalar product (operator *)");

            sparseVectorB.set(0, 3);
            assertEqual(45, vector * sparseVectorB, "Scalar product, modified (operator *)");
        }

};

class VectorNormTest : public GeNuSys::Tests::TestSuite
{

    public:

        VectorNormTest(): TestSuite("VectorNorm") {}

        void run()
        {
            GeNuSys::LinAlg::Vector<int> unit(3);
            unit.set(0, 1);
            unit.set(1, 0);
            unit.set(2, 0);

            std::cout << "Testing unit vector " << unit << std::endl;

            assertEqual(1, GeNuSys::LinAlg::PNorm<1>::norm(unit), "1-norm of unit vector");
            assertEqual(1.0, GeNuSys::LinAlg::PNorm<2>::norm(unit), "2-norm of unit vector");
            assertEqual(1, GeNuSys::LinAlg::PNorm<00>::norm(unit), "infinity norm of unit vector");

            GeNuSys::LinAlg::SparseVector<int> sp_unit(3);
            sp_unit.set(1, 1);

            std::cout << "Testing sparse unit vector " << sp_unit << std::endl;

            assertEqual(1, GeNuSys::LinAlg::PNorm<1>::norm(sp_unit), "1-norm of sparse unit vector");
            assertEqual(1.0, GeNuSys::LinAlg::PNorm<2>::norm(sp_unit), "2-norm of sparse unit vector");
            assertEqual(1, GeNuSys::LinAlg::PNorm<00>::norm(sp_unit), "infinity norm of sparse unit vector");

            GeNuSys::LinAlg::Vector<int> vec1(4);
            vec1.set(0, -2);
            vec1.set(1, 2);
            vec1.set(2, -2);
            vec1.set(3, 2);

            std::cout << "Testing vector " << vec1 << std::endl;

            assertEqual(8, GeNuSys::LinAlg::PNorm<1>::norm(vec1), "1-norm of vector");
            assertEqual(4.0, GeNuSys::LinAlg::PNorm<2>::norm(vec1), "2-norm of vector");
            assertEqual(2, GeNuSys::LinAlg::PNorm<00>::norm(vec1), "infinity norm of vector");

            GeNuSys::LinAlg::SparseVector<int> sp_vec1(8);
            sp_vec1.set(0, -2);
            sp_vec1.set(3, 2);
            sp_vec1.set(5, -2);
            sp_vec1.set(6, 2);

            std::cout << "Testing sparse vector " << sp_vec1 << std::endl;

            assertEqual(8, GeNuSys::LinAlg::PNorm<1>::norm(sp_vec1), "1-norm of sparse vector");
            assertEqual(4.0, GeNuSys::LinAlg::PNorm<2>::norm(sp_vec1), "2-norm of sparse vector");
            assertEqual(2, GeNuSys::LinAlg::PNorm<00>::norm(sp_vec1), "infinity norm of sparse vector");

            GeNuSys::LinAlg::SparseVector<int> empty_sp_vec(16);

            std::cout << "Testing empty sparse vector " << empty_sp_vec << std::endl;

            assertEqual(0, GeNuSys::LinAlg::PNorm<1>::norm(empty_sp_vec), "1-norm of empty sparse vector");
            assertEqual(0.0, GeNuSys::LinAlg::PNorm<2>::norm(empty_sp_vec), "2-norm of empty sparse vector");
            assertEqual(0, GeNuSys::LinAlg::PNorm<00>::norm(empty_sp_vec), "infinity norm of empty sparse vector");
        }

};

class MatrixTest : public GeNuSys::Tests::TestSuite
{

    public:

        MatrixTest(): TestSuite("Matrix") {}

        void run()
        {
            GeNuSys::LinAlg::Matrix<int> matrix(4, 4);

            assertEqual<unsigned int>(4, matrix.getCols(), "Number of matrix cols match");

            assertEqual<unsigned int>(4, matrix.getRows(), "Number of matrix rows match");

            assertEqual(0, matrix(2, 3), "Matrix initialization");

            matrix.set(2, 3, 1);
            matrix.set(0, 2, 1);

            assertEqual(1, matrix(2, 3), "Matrix set and access");
            assertEqual(1, matrix(0, 2), "Matrix set and access II");

            GeNuSys::LinAlg::Matrix<int> copyMatrix = matrix;

            assertEqual(matrix.getCols(), copyMatrix.getCols(), "Number of cols in matrix copy match");

            assertEqual(matrix.getRows(), copyMatrix.getRows(), "Number of rows in matrix copy match");

            bool copyEquals = true;
            for (unsigned int i = 0; i < copyMatrix.getRows() && copyEquals; ++i)
            {
                for (unsigned int j = 0; j < copyMatrix.getCols() && copyEquals; ++j)
                {
                    copyEquals = (matrix(i, j) == copyMatrix(i, j));
                }
            }

            assertTrue(copyEquals, "Matrix copy matches original");

            GeNuSys::LinAlg::Matrix<int> largerMatrix(6, 8);
            largerMatrix.set(4, 1, 1);
            largerMatrix.set(3, 5, 1);
            largerMatrix.set(5, 7, 1);

            matrix = largerMatrix;

            assertEqual(largerMatrix.getCols(), matrix.getCols(), "Number of cols after matrix assign match");

            assertEqual(largerMatrix.getRows(), matrix.getRows(), "Number of rows after matrix assign match");

            bool assignEquals = true;
            for (unsigned int i = 0; i < matrix.getRows() && assignEquals; ++i)
            {
                for (unsigned int j = 0; j < matrix.getCols() && assignEquals; ++j)
                {
                    assignEquals = (matrix(i, j) == largerMatrix(i, j));
                }
            }

            assertTrue(assignEquals, "Matrix matches original after assign");

            GeNuSys::LinAlg::Matrix<double> doubleMatrix = matrix;

            assertEqual(matrix.getCols(), doubleMatrix.getCols(), "Number of cols in converted matrix copy match");

            assertEqual(matrix.getRows(), doubleMatrix.getRows(), "Number of rows in converted matrix copy match");

            bool convertCopyEquals = true;
            for (unsigned int i = 0; i < doubleMatrix.getRows() && convertCopyEquals; ++i)
            {
                for (unsigned int j = 0; j < doubleMatrix.getCols() && convertCopyEquals; ++j)
                {
                    convertCopyEquals = (doubleMatrix(i, j) == (double) matrix(i, j));
                }
            }

            assertTrue(convertCopyEquals, "Converted matrix copy matches original");

            doubleMatrix = largerMatrix;

            assertEqual(largerMatrix.getCols(), doubleMatrix.getCols(), "Number of cols in converted matrix match after assign");

            assertEqual(largerMatrix.getRows(), doubleMatrix.getRows(), "Number of rows in converted matrix match after assign");

            bool convertAssignEquals = true;
            for (unsigned int i = 0; i < doubleMatrix.getRows() && convertAssignEquals; ++i)
            {
                for (unsigned int j = 0; j < doubleMatrix.getCols() && convertAssignEquals; ++j)
                {
                    convertAssignEquals = (doubleMatrix(i, j) == (double) largerMatrix(i, j));
                }
            }

            assertTrue(convertAssignEquals, "Converted matrix matches original after assign");

            GeNuSys::LinAlg::SparseMatrix<int> sparseMatrix(4, 4);
            sparseMatrix.set(1, 2, 1);
            sparseMatrix.set(3, 1, 1);

            GeNuSys::LinAlg::Matrix<int> matrixB = sparseMatrix;

            assertEqual(sparseMatrix.getCols(), matrixB.getCols(), "Number of cols in matrix match after copy sparse");

            assertEqual(sparseMatrix.getRows(), matrixB.getRows(), "Number of rows in matrix match after copy sparse");

            bool copySparseEquals = true;
            for (unsigned int i = 0; i < matrixB.getRows() && copySparseEquals; ++i)
            {
                for (unsigned int j = 0; j < matrixB.getCols() && copySparseEquals; ++j)
                {
                    copySparseEquals = (matrixB(i, j) == sparseMatrix(i, j));
                }
            }

            assertTrue(copySparseEquals, "Matrix matches original after copy sparse");

            GeNuSys::LinAlg::SparseMatrix<int> largerSparseMatrix(10, 7);
            largerSparseMatrix.set(5, 6, 1);
            largerSparseMatrix.set(1, 3, 1);
            largerSparseMatrix.set(4, 5, 1);

            matrixB = largerSparseMatrix;

            assertEqual(largerSparseMatrix.getCols(), matrixB.getCols(), "Number of cols in matrix match after assign sparse");

            assertEqual(largerSparseMatrix.getRows(), matrixB.getRows(), "Number of rows in matrix match after assign sparse");

            bool assignSparseEquals = true;
            for (unsigned int i = 0; i < matrixB.getRows() && assignSparseEquals; ++i)
            {
                for (unsigned int j = 0; j < matrixB.getCols() && assignSparseEquals; ++j)
                {
                    assignSparseEquals = (matrixB(i, j) == largerSparseMatrix(i, j));
                }
            }

            assertTrue(assignSparseEquals, "Matrix matches original after assign sparse");
        }

};

class MatrixNormTest : public GeNuSys::Tests::TestSuite
{

    public:

        MatrixNormTest(): TestSuite("MatrixNorm") {}

        void run()
        {
            GeNuSys::LinAlg::Matrix<int> unit = GeNuSys::LinAlg::Matrix<int>::identity(4, 4);

            std::cout << "Testing identity matrix " << unit << std::endl;

            assertEqual(1, GeNuSys::LinAlg::PNorm<1>::norm(unit), "1-norm of identity matrix");
            assertEqual(1, GeNuSys::LinAlg::PNorm<00>::norm(unit), "infinity norm of identity matrix");
            assertEqual(2.0, GeNuSys::LinAlg::FrobeniusNorm::norm(unit), "Frobenius norm of identity matrix");

            GeNuSys::LinAlg::SparseMatrix<int> sp_unit = GeNuSys::LinAlg::SparseMatrix<int>::identity(4, 4);

            std::cout << "Testing sparse unit vector " << sp_unit << std::endl;

            assertEqual(1, GeNuSys::LinAlg::PNorm<1>::norm(sp_unit), "1-norm of sparse identity matrix");
            assertEqual(1, GeNuSys::LinAlg::PNorm<00>::norm(sp_unit), "infinity norm of sparse identity matrix");
            assertEqual(2.0, GeNuSys::LinAlg::FrobeniusNorm::norm(sp_unit), "Frobenius norm of sparse identity matrix");

            GeNuSys::LinAlg::Matrix<int> zeros(4, 4);

            std::cout << "Testing empty matrix " << zeros << std::endl;

            assertEqual(0, GeNuSys::LinAlg::PNorm<1>::norm(zeros), "1-norm of empty matrix");
            assertEqual(0, GeNuSys::LinAlg::PNorm<00>::norm(zeros), "infinity norm of empty matrix");
            assertEqual(0.0, GeNuSys::LinAlg::FrobeniusNorm::norm(zeros), "Frobenius norm of empty matrix");

            GeNuSys::LinAlg::SparseMatrix<int> sp_zeros(4, 4);

            std::cout << "Testing empty sparse matrix " << sp_zeros << std::endl;

            assertEqual(0, GeNuSys::LinAlg::PNorm<1>::norm(sp_zeros), "1-norm of empty sparse matrix");
            assertEqual(0, GeNuSys::LinAlg::PNorm<00>::norm(sp_zeros), "infinity norm of empty sparse matrix");
            assertEqual(0.0, GeNuSys::LinAlg::FrobeniusNorm::norm(sp_zeros), "Frobenius norm of empty sparse matrix");

            GeNuSys::LinAlg::Matrix<int> mat1(4, 4);
            for (unsigned int i = 0; i < 4; ++i)
            {
                for (unsigned int j = 0; j < 4; ++j)
                {
                    if (i <= j)
                    {
                        mat1.set(i, j, -((int)i + 2));
                    }
                    else
                    {
                        mat1.set(i, j, j + 1);
                    }
                }
            }

            std::cout << "Testing matrix " << mat1 << std::endl;

            assertEqual(14, GeNuSys::LinAlg::PNorm<1>::norm(mat1), "1-norm of matrix");
            assertEqual(11, GeNuSys::LinAlg::PNorm<00>::norm(mat1), "infinity norm of matrix");
            double fnorm_mat1 = GeNuSys::LinAlg::FrobeniusNorm::norm(mat1);
            assertDifference(120.0, 0.0001, fnorm_mat1 * fnorm_mat1, "Frobenius norm of matrix");

            GeNuSys::LinAlg::SparseMatrix<int> sp_mat1(4, 4);
            sp_mat1.set(0, 0, 1);
            sp_mat1.set(1, 0, -1);
            sp_mat1.set(1, 2, -2);
            sp_mat1.set(2, 0, 2);
            sp_mat1.set(3, 1, 3);
            sp_mat1.set(2, 3, -4);
            sp_mat1.set(3, 3, 4);

            std::cout << "Testing sparse matrix " << sp_mat1 << std::endl;

            assertEqual(8, GeNuSys::LinAlg::PNorm<1>::norm(sp_mat1), "1-norm of sparse matrix");
            assertEqual(7, GeNuSys::LinAlg::PNorm<00>::norm(sp_mat1), "infinity norm of sparse matrix");
            double fnorm_sp_mat1 = GeNuSys::LinAlg::FrobeniusNorm::norm(sp_mat1);
            assertDifference(51.0, 0.0001, fnorm_sp_mat1 * fnorm_sp_mat1, "Frobenius norm of sparse matrix");
        }

};

/* Hátralévő
 * matrix::diag(vec)
 *    template<typename ElementType>
    static void vct_mod(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result);

    template<typename ElementType>
    static void vct_mod(Vector<ElementType>& vct, const ElementType& value);

    template<typename ElementType>
    static void vct_mod(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result);

    template<typename ElementType>
    static void vct_mod(SparseVector<ElementType>& vct, const ElementType& value);

    template<typename ElementType>
    static void vct_mods(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result);

    template<typename ElementType>
    static void vct_mods(Vector<ElementType>& vct, const ElementType& value);

    template<typename ElementType>
    static void vct_mods(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result);

    template<typename ElementType>
    static void vct_mods(SparseVector<ElementType>& vct, const ElementType& value);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(Matrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(SparseMatrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_idiv(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_idiv(Matrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_idiv(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_idiv(SparseMatrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_mod(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mod(Matrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_mod(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mod(SparseMatrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_mods(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mods(Matrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_mods(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mods(SparseMatrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_div(const Matrix<ElementType>& mat, const ElementType& value, Matrix<typename ElementTraits<ElementType>::RationalType>& result);

    template<typename ElementType>
    static void mat_div(Matrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_div(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<typename ElementTraits<ElementType>::RationalType>& result);

    template<typename ElementType>
    static void mat_div(SparseMatrix<ElementType>& mat, const ElementType& value);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& mat, const Vector<ElementType>& vct, Vector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& mat, const Vector<ElementType>& vct, SparseVector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& mat, const SparseVector<ElementType>& vct, Vector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& mat, const SparseVector<ElementType>& vct, SparseVector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& mat, const Vector<ElementType>& vct, Vector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& mat, const Vector<ElementType>& vct, SparseVector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& mat, const SparseVector<ElementType>& vct, Vector<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& mat, const SparseVector<ElementType>& vct, SparseVector<ElementType>& result);

    template<typename ElementType>
    static void mat_add(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_add(Matrix<ElementType>& op1, const Matrix<ElementType>& op2);

    template<typename ElementType>
    static void mat_add(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_add(Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

    template<typename ElementType>
    static void mat_add(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_add(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_sub(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_sub(Matrix<ElementType>& op1, const Matrix<ElementType>& op2);

    template<typename ElementType>
    static void mat_sub(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_sub(Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

    template<typename ElementType>
    static void mat_sub(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_sub(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

    template<typename ElementType>
    static void mat_mul(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

    template<typename ElementType>
    static bool mat_eq(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2);

    template<typename ElementType>
    static bool mat_eq(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

    template<typename ElementType>
    static bool mat_eq(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2);

    template<typename ElementType>
    static bool mat_eq(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);
//    getsubvector,getrow,swapcols,swaprows,getcol,getsubmatrix (linalg_traits.h)
 */
/*
        static ElementType sgn(const ElementType& value);

        static typename ElementTypeTraits<ElementType>::AbsType abs(const ElementType& value);

        static typename ElementTypeTraits<ElementType>::AbsSqrType absSqr(const ElementType& value);

        static ElementType conj(const ElementType& value);

        static ElementType pow(const ElementType& value, unsigned int n);

        static typename ElementTypeTraits<ElementType>::RealType sqrt(const ElementType& value);

        static typename ElementTypeTraits<ElementType>::RealType root(const ElementType& value, int n);

        static typename ElementTypeTraits<ElementType>::RationalType div(const ElementType& a, const ElementType& b);

        static bool divisible(const ElementType& a, const ElementType& b);

        static ElementType idiv(const ElementType& a, const ElementType& b);

        static ElementType mod(const ElementType& a, const ElementType& b);

        static ElementType mods(const ElementType& a, const ElementType& b);
 */
