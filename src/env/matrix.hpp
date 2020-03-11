/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2018 by David B. Blumenthal                              *
 *                                                                          *
 *   This file is part of GEDLIB.                                           *
 *                                                                          *
 *   GEDLIB is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 3 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   GEDLIB is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                          *
 ***************************************************************************/

/*!
 * @file matrix.hpp
 * @brief ged::Matrix class declaration.
 */

#ifndef SRC_ENV_MATRIX_HPP_
#define SRC_ENV_MATRIX_HPP_

#include "common_types.hpp"

namespace ged {

/*!
 * @brief A matrix class with basic functionality.
 */
template<class ScalarT>
class Matrix {

public:
	/*!
	 * @brief Constructs an empty cost matrix.
	 */
	Matrix();

	/*!
	 * @brief Constructs a matrix of given size with default entries.
	 * @param[in] num_rows The number of rows.
	 * @param[in] num_cols The number of columns.
	 * @param[in] val The default value of all cells.
	 */
	Matrix(std::size_t num_rows, std::size_t num_cols, ScalarT val = 0);

	/*!
	 * @brief Copy constructor.
	 * @param[in] matrix Matrix that should be copied.
	 */
	Matrix(const Matrix & matrix);

	/*!
	 * @brief Assignment operator.
	 * @param[in] matrix Matrix that should be assigned.
	 */
	void operator=(const Matrix & matrix);

	/*!
	 * @brief Swaps the matrix with another matrix.
	 * @param[in,out] matrix Matrix that is to be swapped with the calling matrix.
	 */
	void swap(Matrix& matrix);

	/*!
	 * @brief Provides access to a cell.
	 * @param[in] row Row of the cell.
	 * @param[in] col Column of the cell.
	 * @return Constant reference to the cell (\p row,\p col).
	 */
	const ScalarT & operator() (std::size_t row, std::size_t col) const;

	/*!
	 * @brief Provides access to a cell.
	 * @param[in] row Row of the cell.
	 * @param[in] col Column of the cell.
	 * @return Reference to the cell (\p row,\p col).
	 */
	ScalarT & operator() (std::size_t row, std::size_t col);


	/*!
	 * @brief Provides access to internal data.
	 * @return Pointer to the data array.
	 */
	ScalarT * data();

	/*!
	 * @brief Provides constant access to internal data.
	 * @return Constant pointer to the data array.
	 */
	const ScalarT * data() const;


	/*!
	 * @brief Returns the number of rows.
	 * @return The number of rows.
	 */
	std::size_t num_rows() const;

	/*!
	 * @brief Returns the number of columns.
	 * @return The number of columns.
	 */
	std::size_t num_cols() const;

	/*!
	 * @brief Resizes the matrix.
	 * @param[in] num_rows New number of rows.
	 * @param[in] num_cols New number of columns.
	 */
	void resize(std::size_t num_rows, std::size_t num_cols);

	/*!
	 * \brief Sets all cells to \p val.
	 */
	void set_to_val(const ScalarT & val);

	/*!
	 * \brief Takes the matrix to the power of \p n.
	 */
	void power(std::size_t n);

	/*!
	 * \brief Returns the maximal coefficient.
	 */
	ScalarT max() const;

	/*!
	 * \brief Returns the minimal coefficient.
	 */
	ScalarT min() const;

	/*!
	 * \brief Transposes the matrix.
	 */
	void transpose();

	/*!
	 * \brief Returns the transposed matrix.
	 */
	Matrix<ScalarT> transposed() const;

	/*!
	 * \brief Matrix-scalar multiplication assignment.
	 */
	Matrix<ScalarT> & operator*=(const ScalarT & scalar);

	/*!
	 * \brief Matrix-scalar multiplication assignment.
	 */
	Matrix<ScalarT> & operator/=(const ScalarT & scalar);

	/*!
	 * \brief Matrix-matrix addition assignment.
	 */
	Matrix<ScalarT> & operator+=(const Matrix<ScalarT> & matrix);

	/*!
	 * \brief Matrix-matrix substraction assignment.
	 */
	Matrix<ScalarT> & operator-=(const Matrix<ScalarT> & matrix);

	/*!
	 * \brief Matrix-scalar multiplication.
	 */
	Matrix<ScalarT> operator*(const ScalarT & scalar) const;

	/*!
	 * \brief Matrix-scalar multiplication.
	 */
	Matrix<ScalarT> operator/(const ScalarT & scalar) const;

	/*!
	 * \brief Matrix-matrix addition.
	 */
	Matrix<ScalarT> operator+(const Matrix<ScalarT> & matrix) const;

	/*!
	 * \brief Matrix-matrix addition.
	 */
	Matrix<ScalarT> operator-(const Matrix<ScalarT> & matrix) const;

	/*!
	 * \brief Returns reference to the internal Eigen matrix.
	 */
	Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> & matrix();

	/*!
	 * \brief Returns constant reference to the internal Eigen matrix.
	 */
	const Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> & matrix() const;

private:
	Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> matrix_;

};

/*!
 * @brief Matrix with double entries.
 */
typedef Matrix<double> DMatrix;

/*!
 * @brief Matrix with int entries.
 */
typedef Matrix<int> IMatrix;

/*!
 * @brief Streams a matrix.
 * @param[in,out] os Output stream.
 * @param[in] matrix The matrix that should be streamed.
 * @return The output stream @p os.
 * @relates ged::Matrix
 */
template<typename ScalarT>
std::ostream & operator<<(std::ostream & os, const Matrix<ScalarT> & matrix);

}

#include "matrix.ipp"

#endif /* SRC_ENV_MATRIX_HPP_ */
