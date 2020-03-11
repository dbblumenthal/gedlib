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
 * \file matrix.ipp
 * \brief ged::Matrix class definition.
 */

#ifndef MATRIX_IPP
#define MATRIX_IPP 1

namespace ged {

template<class ScalarT>
Matrix<ScalarT>::
Matrix(std::size_t num_rows, std::size_t num_cols, ScalarT val) :
matrix_(Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>::Constant(num_rows, num_cols, val)) {}

template<class ScalarT>
Matrix<ScalarT>::
Matrix() :
matrix_(){}

template<class ScalarT>
Matrix<ScalarT>::
Matrix(const Matrix<ScalarT> & matrix) :
matrix_(matrix.matrix_) {}

template<class ScalarT>
void
Matrix<ScalarT>::
operator=(const Matrix & matrix) {
	matrix_ = matrix.matrix_;
}

template<class ScalarT>
const ScalarT &
Matrix<ScalarT>::
operator() (std::size_t row, std::size_t col) const {
	return matrix_(row, col);
}

template<class ScalarT>
ScalarT &
Matrix<ScalarT>::
operator() (std::size_t row, std::size_t col) {
	return matrix_(row, col);
}

template<class ScalarT>
ScalarT *
Matrix<ScalarT>::
data() {
	return matrix_.data();
}

template<class ScalarT>
const ScalarT *
Matrix<ScalarT>::
data() const {
	return matrix_.data();
}

template<class ScalarT>
std::size_t
Matrix<ScalarT>::
num_rows() const {
	return static_cast<std::size_t>(matrix_.rows());
}

template<class ScalarT>
void
Matrix<ScalarT>::
resize(std::size_t num_rows, std::size_t num_cols) {
	matrix_.resize(num_rows, num_cols);
}

template<class ScalarT>
void
Matrix<ScalarT>::
set_to_val(const ScalarT & val) {
	for (std::size_t row{0}; row < num_rows(); row++) {
		for (std::size_t col{0}; col < num_cols(); col++) {
			(*this)(row, col) = val;
		}
	}
}

template<class ScalarT>
std::size_t
Matrix<ScalarT>::
num_cols() const {
	return static_cast<std::size_t>(matrix_.cols());
}

template<class ScalarT>
void
Matrix<ScalarT>::
power(std::size_t n) {
	if (n <= 1) {
		return;
	}
	Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> temp = Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>::Ones(num_rows(), num_cols());
	while (n > 1) {
		if (n % 2 == 0) {
			matrix_ *= matrix_;
			n /= 2;
		}
		else {
			temp *= matrix_;
			matrix_ *= matrix_;
			n--;
			n /= 2;
		}
	}
	matrix_ *= temp;
}

template<class ScalarT>
void
Matrix<ScalarT>::
transpose() {
	matrix_.transpose();
}

template<class ScalarT>
Matrix<ScalarT>
Matrix<ScalarT>::
transposed() const {
	Matrix<ScalarT> transposed_matrix(*this);
	transposed_matrix.transpose();
	return transposed_matrix;
}

template<class ScalarT>
void
Matrix<ScalarT>::
swap(Matrix & rhs) {
	std::swap(matrix_, rhs.matrix_);
}

template<class ScalarT>
ScalarT
Matrix<ScalarT>::
max() const {
	return matrix_.maxCoeff();
}

template<class ScalarT>
ScalarT
Matrix<ScalarT>::
min() const {
	return matrix_.minCoeff();
}

template<class ScalarT>
Matrix<ScalarT> &
Matrix<ScalarT>::
operator*=(const ScalarT & scalar) {
	matrix_ *= scalar;
	return (*this);
}

template<class ScalarT>
Matrix<ScalarT> &
Matrix<ScalarT>::
operator/=(const ScalarT & scalar) {
	matrix_ /= scalar;
	return (*this);
}

template<class ScalarT>
Matrix<ScalarT> &
Matrix<ScalarT>::
operator+=(const Matrix<ScalarT> & matrix) {
	matrix_ += matrix.matrix_;
	return (*this);
}

template<class ScalarT>
Matrix<ScalarT> &
Matrix<ScalarT>::
operator-=(const Matrix<ScalarT> & matrix) {
	matrix_ -= matrix.matrix_;
	return (*this);
}

template<class ScalarT>
Matrix<ScalarT>
Matrix<ScalarT>::
operator*(const ScalarT & scalar) const {
	Matrix<ScalarT> new_matrix(*this);
	new_matrix.matrix_ *= scalar;
	return new_matrix;
}

template<class ScalarT>
Matrix<ScalarT>
Matrix<ScalarT>::
operator/(const ScalarT & scalar) const {
	Matrix<ScalarT> new_matrix(*this);
	new_matrix.matrix_ /= scalar;
	return new_matrix;
}

template<class ScalarT>
Matrix<ScalarT>
Matrix<ScalarT>::
operator+(const Matrix<ScalarT> & matrix) const {
	Matrix<ScalarT> new_matrix(*this);
	new_matrix.matrix_ += matrix.matrix_;
	return new_matrix;
}

template<class ScalarT>
Matrix<ScalarT>
Matrix<ScalarT>::
operator-(const Matrix<ScalarT> & matrix) const {
	Matrix<ScalarT> new_matrix(*this);
	new_matrix.matrix_ -= matrix.matrix_;
	return new_matrix;
}

template<class ScalarT>
Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> &
Matrix<ScalarT>::
matrix() {
	return matrix_;
}

template<class ScalarT>
const Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> &
Matrix<ScalarT>::
matrix() const {
	return matrix_;
}

template<class ScalarT>
std::ostream & operator<<(std::ostream & os, const Matrix<ScalarT> & matrix) {
	os << matrix.matrix();
	return os;
}

}

#endif
