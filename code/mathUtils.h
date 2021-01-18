#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <type_traits>
#include <mkl.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef int32_t s32;
typedef int16_t s16;
typedef double r64;
typedef uint32_t b32;

#define PI64 3.14159265358979323846 
#define IPI64 0.318309886183790671538 
#define IS_EVEN(n) (((n) % 2) == 0)
#define MIN(a,b) ((a) > (b) ? (b) : (a))

#define count(array) (sizeof(array)/sizeof((array)[0]))

#ifdef CI_DEBUG
	#define assert(expr) if(!(expr)){*(int*)0 = 0;}
#else
	#define assert(expr)
#endif

struct Matrix;
struct Vector;

enum MatrixType
{
	MATRIXTYPE_GENERAL,
	MATRIXTYPE_VECTOR,
	MATRIXTYPE_IMPLICIT,
};

enum OperationType
{
	OPERATIONTYPE_SUM,
	OPERATIONTYPE_SUBSTRACTION,
	OPERATIONTYPE_PRODUCT,
	OPERATIONTYPE_IMPLICIT,
};

enum ExpressionState
{
	EXPRESSIONSTATE_UNINITIALIZED,
	EXPRESSIONSTATE_INITIALIZED,
};

enum MATRIX_OP_FLAG
{
	M_ROW_MAJOR = 0,
	M_COLUMN_MAJOR = 1,
	M_NO_TRANSPOSE = 2,
	M_TRANSPOSE = 4,
	M_ADJOINT = 8,
	M_CONJUGATE = 16,

	M_NO_FLAG = 64
};
int cblasFlagTable[] = {CblasRowMajor, CblasColMajor, CblasNoTrans, 
								M_NO_FLAG, CblasTrans, M_NO_FLAG, 
								M_NO_FLAG, M_NO_FLAG, CblasConjTrans};

char mklFlagTable[] = {'R', 'C', 'N', M_NO_FLAG, 'T', M_NO_FLAG, 
							  M_NO_FLAG, M_NO_FLAG, 'C', M_NO_FLAG, 
							  M_NO_FLAG, M_NO_FLAG, M_NO_FLAG, M_NO_FLAG, 
							  M_NO_FLAG, M_NO_FLAG, 'R'};

inline void scaleAndAddMatrices(Matrix& result, Matrix m1, Matrix m2, r64 m1Coefficient, r64 m2Coefficient, 
										  MATRIX_OP_FLAG m1Transpose, MATRIX_OP_FLAG m2Transpose);

inline void scaleAndMultiplyMatrices(Matrix& result, Matrix m1, Matrix m2, r64 m1m2Coefficient, r64 resultCoefficient, 
                				  			 MATRIX_OP_FLAG m1Transpose, MATRIX_OP_FLAG m2Transpose);

inline void scaleAndMultiplyMatrices(Matrix& result, Matrix m1, Matrix m2, r64 m1m2Coefficient, 
												 MATRIX_OP_FLAG m1Transpose, MATRIX_OP_FLAG m2Transpose);

//inline Vector vmmul(Vector v, Matrix m);
//inline void madd(Matrix& result, Matrix m1, Matrix m2);
//inline void mmuladd(Matrix& result, Matrix m1, Matrix m2);
//inline void mmul(Matrix& result, Matrix m1, Matrix m2, bool, bool);
//inline Matrix mmul(Matrix m1, Matrix m2, bool, bool);
//inline void mmul(Matrix& result, Matrix m1, Matrix m2);
//inline Matrix mmul(Matrix m1, Matrix m2);
inline int orthonormalizeVectors(Matrix basis, r64* workArea);

template<typename A, typename B>
struct IsSame
{
	enum {value = 0};
};

template<typename A>
struct IsSame<A,A>
{
	enum {value = 1};
};

template <typename T>
struct remove_cvref 
{
    typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

#define DEFINE_HAS_MEMBER_FOR_ID(name) \
template <typename T, typename = int> \
struct Has_##name { enum {value = 0}; }; \
template <typename T> \
struct Has_##name <T, decltype((void)T::##name, 0)> { enum {value = 1}; };

DEFINE_HAS_MEMBER_FOR_ID(leftExpr)
DEFINE_HAS_MEMBER_FOR_ID(rightExpr)

template <typename Expression>
constexpr bool isExpandableExpression = (Has_leftExpr<Expression>::value ||
													  Has_rightExpr<Expression>::value);

struct MatrixOperationMarker {};


struct transposed : public MatrixOperationMarker
{
	const Matrix& m;

	transposed(const Matrix& matrix) : m(matrix){}
};

template <typename Expression1, typename Expression2, OperationType /*opType*/type>
struct MatrixOperation : public MatrixOperationMarker
{
	const Expression1& leftExpr;	
	const Expression2& rightExpr;
	enum {opType = type};
	const r64 leftCoefficient;
	const r64 rightCoefficient;
	MATRIX_OP_FLAG leftFlag;
	MATRIX_OP_FLAG rightFlag;

	MatrixOperation(const Expression1& left, const Expression2& right) : leftExpr(left), 
						 rightExpr(right), leftFlag(M_NO_TRANSPOSE), rightFlag(M_NO_TRANSPOSE), leftCoefficient(1), rightCoefficient(1){}
	MatrixOperation(const Expression1& left, const Expression2& right, MATRIX_OP_FLAG lFlag, 
						 MATRIX_OP_FLAG rFlag, r64 lSign = 1, r64 rSign = 1) : leftExpr(left), rightExpr(right), leftFlag(lFlag), rightFlag(rFlag), leftCoefficient(lSign), rightCoefficient(rSign){}
};


/* ### Type restriction of suitable operands in math expressions inspired by https://www.youtube.com/watch?v=4IUCBx5fIv0 ### */
template <typename Op>
constexpr bool isSuitableOperand = ((IsSame<Op, Matrix>::value == 1) || 
												(IsSame<Op, Vector>::value == 1) || 
												std::is_base_of_v<MatrixOperationMarker, std::remove_reference_t<Op>>); 

template <typename LHS, typename RHS>
constexpr bool areSuitableOperands = ((isSuitableOperand<LHS> && isSuitableOperand<RHS>) || 
												  (std::is_scalar_v<LHS> && isSuitableOperand<RHS>) ||
												  (std::is_scalar_v<RHS> && isSuitableOperand<LHS>)) &&
												  (!IsSame<LHS, Vector>::value || !IsSame<RHS, Vector>::value);

template <typename T, int used>
struct SeparatedScalarProduct
{
	const r64 scalar;
	const T& other;
	enum {isUsed = used};
	SeparatedScalarProduct(r64 s, const T& multiplicant) : scalar(s), other(multiplicant) {}
};

template <typename Expression1, typename Expression2, OperationType type>
constexpr auto separateScalarCoefficient(const MatrixOperation<Expression1, Expression2, type>& op)
{
	SeparatedScalarProduct<MatrixOperation<Expression1, Expression2, type>, 0> result(1.0, op);
	return result;
}

template <typename Expression1, typename Expression2, typename = std::enable_if_t<std::is_scalar_v<Expression1>>, typename = std::enable_if_t<!std::is_scalar_v<Expression2>>>
constexpr auto separateScalarCoefficient(const MatrixOperation<Expression1, Expression2, OPERATIONTYPE_PRODUCT>& op)
{
	SeparatedScalarProduct<Expression2, 1> result((r64)op.leftExpr, op.rightExpr);
	return result;
}

template <typename Expression1, typename Expression2, typename = std::enable_if_t<std::is_scalar_v<Expression2>>>
constexpr auto separateScalarCoefficient(const MatrixOperation<Expression1, Expression2, OPERATIONTYPE_PRODUCT>& op)
{
	SeparatedScalarProduct<Expression1, 1> result((r64)op.rightExpr, op.leftExpr);
	return result;
}

template <typename Expression1, typename Expression2, OperationType type>
constexpr auto stripTraitAndConstructOperation(const Expression1& e1, const Expression2& e2)
{
	if constexpr(isExpandableExpression<Expression1>)
	{
		const SeparatedScalarProduct e1Separated = separateScalarCoefficient(e1);

		if constexpr(isExpandableExpression<Expression2>)
		{
			const SeparatedScalarProduct e2Separated = separateScalarCoefficient(e2);

			MatrixOperation<remove_cvref_t<decltype(e1Separated.other)>, remove_cvref_t<decltype(e2Separated.other)>, type> 
				result(e1Separated.other, e2Separated.other, M_NO_TRANSPOSE, M_NO_TRANSPOSE, e1Separated.scalar, e2Separated.scalar);
			return result;
		}
		else
		{
			MatrixOperation<remove_cvref_t<decltype(e1Separated.other)>, Expression2, type> 
				result(e1Separated.other, e2, M_NO_TRANSPOSE, M_NO_TRANSPOSE, e1Separated.scalar, 1.0);
			return result;
		}
	} 
	else if constexpr(IsSame<Expression1, transposed>::value)
	{
		if constexpr(IsSame<Expression2, transposed>::value)
		{
			MatrixOperation<Matrix, Matrix, type> result(e1.m, e2.m, M_TRANSPOSE, M_TRANSPOSE);
			return result;
		}

		MatrixOperation<Matrix, Expression2, type> result(e1.m, e2, M_TRANSPOSE, M_NO_TRANSPOSE);
		return result;
	}
	else if constexpr(IsSame<Expression2, transposed>::value)
	{
		MatrixOperation<Expression1, Matrix, type> result(e1, e2.m, M_NO_TRANSPOSE, M_TRANSPOSE);
		return result;
	}
	else
	{
		MatrixOperation<Expression1, Expression2, type> result(e1, e2);
		return result;
	}

}

template <typename Expression1, typename Expression2, typename = std::enable_if_t<areSuitableOperands<Expression1, Expression2>>> 
inline auto operator+(const Expression1& e1, const Expression2& e2)
{
	auto result = stripTraitAndConstructOperation<Expression1, Expression2, OPERATIONTYPE_SUM>(e1, e2);

	return result;
}


template <typename Expression1, typename Expression2, typename = std::enable_if_t<areSuitableOperands<Expression1, Expression2>>> 
inline auto operator-(const Expression1& e1, const Expression2& e2)
{
	auto collapsedResult = stripTraitAndConstructOperation<Expression1, Expression2, OPERATIONTYPE_SUM>(e1, e2);

	MatrixOperation<remove_cvref_t<decltype(collapsedResult.leftExpr)>, remove_cvref_t<decltype(collapsedResult.rightExpr)>, OPERATIONTYPE_SUM> 
		result(collapsedResult.leftExpr, collapsedResult.rightExpr, M_NO_TRANSPOSE, M_NO_TRANSPOSE, collapsedResult.leftCoefficient, -collapsedResult.rightCoefficient);

	return result;
}

template <typename Expression, typename = std::enable_if_t<isSuitableOperand<Expression>>> 
inline auto operator-(const Expression& e)
{
	auto result = stripTraitAndConstructOperation<r64, Expression, OPERATIONTYPE_PRODUCT>(-1.0, e);

	return result;
}

template <typename Expression1, typename Expression2, typename = std::enable_if_t<areSuitableOperands<Expression1, Expression2>>>
inline auto operator*(const Expression1& e1, const Expression2& e2)
{

	if constexpr(isExpandableExpression<Expression1>)
	{
		static_assert(Expression1::opType != OPERATIONTYPE_PRODUCT,
						  "Multiplication of more than 2 operands(scalars, matrices or vectors) is not supported!");
	}
	if constexpr(isExpandableExpression<Expression2>)
	{
		static_assert(Expression2::opType != OPERATIONTYPE_PRODUCT, 
						  "Multiplication of more than 2 operands(scalars, matrices or vectors) is not supported!");
	}
	
	auto result = stripTraitAndConstructOperation<Expression1, Expression2, OPERATIONTYPE_PRODUCT>(e1, e2);

	return result;
}

template <ExpressionState state, typename Expression1, typename Expression2, OperationType type, typename DestType>
constexpr void expandExpression(const MatrixOperation<Expression1, Expression2, type>& op, DestType& dest)
{
	if constexpr(!isExpandableExpression<Expression1> && !isExpandableExpression<Expression2>)
	{
		expandExpression<state>(op, dest);
	} 

	if constexpr(isExpandableExpression<Expression1>) 
	{
		expandExpression<(!isExpandableExpression<Expression1> && !isExpandableExpression<Expression2>) ? EXPRESSIONSTATE_INITIALIZED : EXPRESSIONSTATE_UNINITIALIZED>
			(op.leftExpr, dest);

		if constexpr(isExpandableExpression<Expression2>)
		{
			expandExpression<EXPRESSIONSTATE_INITIALIZED>(op.rightExpr, dest);
		}
		else
		{
			MatrixOperation<DestType, remove_cvref_t<decltype(op.rightExpr)>, type> result(dest, op.rightExpr);
			expandExpression<EXPRESSIONSTATE_UNINITIALIZED>(result, dest);
		}
	}
	else if constexpr(isExpandableExpression<Expression2>)
	{
		expandExpression<(isExpandableExpression<Expression1>) ? EXPRESSIONSTATE_INITIALIZED : EXPRESSIONSTATE_UNINITIALIZED>
			(op.rightExpr, dest);

		if constexpr(isExpandableExpression<Expression1>)
		{
			expandExpression<EXPRESSIONSTATE_INITIALIZED>(op.leftExpr, dest);
		}
		else
		{
			MatrixOperation<remove_cvref_t<decltype(op.leftExpr)>, DestType, type> result(op.leftExpr, dest);
			expandExpression<EXPRESSIONSTATE_UNINITIALIZED>(result, dest);
		}

	}

}


struct Dim
{
	size_t rows;			
	size_t columns;
};

struct MatrixConversion
{
	Matrix* matrix;
	size_t index;
	MatrixConversion(Matrix* m, size_t i): matrix(m), index(i) {}

	r64& operator[] (size_t row);
	operator r64*();
	operator Vector();
};

struct Matrix
{
	r64* elements;
	union
	{
		Dim dim;
		struct
		{
			size_t rows;			
			size_t columns;
		};
	};

	Matrix(){}

	Matrix(size_t activeRows, size_t activeColumns, r64* memory) : rows(activeRows), columns(activeColumns), elements(memory){}

	Matrix(size_t maxRows, size_t maxColumns) : rows(maxRows), columns(maxColumns)
	{
		elements = (r64*)calloc(maxRows*maxColumns, sizeof(r64));
	}

	Matrix(size_t maxRows, size_t maxColumns, size_t activeRows, size_t activeColumns) : 
		rows(activeRows), columns(activeColumns)
	{
		elements = (r64*)calloc(maxRows*maxColumns, sizeof(r64));
	}

	inline r64& operator[](size_t index)
	{
		return elements[index];
	}

	inline MatrixConversion operator()(size_t row)
	{
		MatrixConversion conversion(this, row);
		return conversion;
	}

	inline r64& operator()(size_t row, size_t column)
	{
		return elements[row*columns + column];
	}

	template <typename Expression1, typename Expression2, OperationType type>
	inline Matrix& operator=(const MatrixOperation<Expression1, Expression2, type>& op)
	{
		expandExpression<EXPRESSIONSTATE_UNINITIALIZED>(op, *this);

		return *this;
	}
	template <typename Expression1, typename Expression2, OperationType type>
	inline Matrix& operator+=(const MatrixOperation<Expression1, Expression2, type>& op)
	{
		expandExpression<EXPRESSIONSTATE_INITIALIZED>(op, *this);

		return *this;
	}

};

struct Vector
{
	r64* elements;
	size_t dim;

	Vector(){}

	//Shallow copy
	Vector(MatrixConversion mc)
	{
		elements = mc.matrix->elements + mc.index*mc.matrix->columns; //A row of matrix
		dim = mc.matrix->columns;
	}

	Vector(size_t size) : dim(size)
	{
		elements = (r64*)calloc(size, sizeof(r64));
	}

	Vector(size_t size, r64* memory) : dim(size), elements(memory){}

	inline r64& operator[](size_t index)
	{
		return elements[index];
	}

	//Deep copy
	inline Vector& operator=(Vector v)
	{
		assert(dim >= v.dim);
		cblas_dcopy(dim, v.elements, 1, elements, 1);
		dim = v.dim;

		return *this;
	}
	
	inline r64 operator*(Vector v)
	{
		assert(dim == v.dim);
		return cblas_ddot(dim, elements, 1, v.elements, 1);
	}

	template <typename Expression1, typename Expression2, OperationType type>
	inline Vector& operator=(const MatrixOperation<Expression1, Expression2, type>& op)
	{
		expandExpression<EXPRESSIONSTATE_UNINITIALIZED>(op, *this);

		return *this;
	}
	template <typename Expression1, typename Expression2, OperationType type>
	inline Vector& operator+=(const MatrixOperation<Expression1, Expression2, type>& op)
	{
		expandExpression<EXPRESSIONSTATE_INITIALIZED>(op, *this);

		return *this;
	}

	inline Vector& operator+=(const Vector& v)
	{
		assert(elements != 0);
		assert(dim == v.dim);

		cblas_daxpy(dim, 1.0, v.elements, 1, elements, 1);

		return *this;
	}
};


// ### Matrix operations ###

//Matrix-Matrix sum, assign the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_UNINITIALIZED>
(const MatrixOperation<Matrix, Matrix, OPERATIONTYPE_SUM>& op, Matrix& dest)
{
	assert(dest.elements != 0);
	scaleAndAddMatrices(dest, op.leftExpr, op.rightExpr, op.leftCoefficient, op.rightCoefficient, op.leftFlag, op.rightFlag);
}

//Matrix-Matrix sum, add the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_INITIALIZED>
(const MatrixOperation<Matrix, Matrix, OPERATIONTYPE_SUM>& op, Matrix& dest)
{
	assert(dest.elements != 0);
	scaleAndAddMatrices(dest, dest, op.leftExpr, 1.0, op.leftCoefficient, M_NO_TRANSPOSE, op.leftFlag); 
	scaleAndAddMatrices(dest, dest, op.rightExpr, 1.0, op.rightCoefficient, M_NO_TRANSPOSE, op.rightFlag);
}

//Matrix-Matrix product, assign the result to this matrix
template<>
inline void expandExpression<EXPRESSIONSTATE_UNINITIALIZED>
(const MatrixOperation<Matrix, Matrix, OPERATIONTYPE_PRODUCT>& op, Matrix& dest)
{
	assert(dest.elements != 0);
	scaleAndMultiplyMatrices(dest, op.leftExpr, op.rightExpr, 1.0, op.leftFlag, op.rightFlag);
}

//Matrix-Matrix product, add the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_INITIALIZED>
(const MatrixOperation<Matrix, Matrix, OPERATIONTYPE_PRODUCT>& op, Matrix& dest)
{
	assert(dest.elements != 0);
	scaleAndMultiplyMatrices(dest, op.leftExpr, op.rightExpr, op.leftCoefficient*op.rightCoefficient, 1.0, op.leftFlag, op.rightFlag);
}

//Scalar-Matrix multiplication, assign or add the result to dest
template <ExpressionState state, typename Expression1, typename = std::enable_if_t<std::is_scalar_v<Expression1>>>
inline void expandExpression(const MatrixOperation<Expression1, Matrix, OPERATIONTYPE_PRODUCT>& op, Matrix& dest)
{
	assert(dest.elements != 0);
	if constexpr(state == EXPRESSIONSTATE_UNINITIALIZED)
	{
		scaleAndAddMatrices(dest, op.rightExpr, op.leftExpr, op.rightFlag);
	}
	else
	{
		scaleAndMultiplyMatrices(dest, dest, op.rightExpr, 1.0, op.leftExpr, M_NO_TRANSPOSE, op.rightFlag);
	}
}
template <ExpressionState state, typename Expression2, typename = std::enable_if_t<std::is_scalar_v<Expression2>>>
inline void expandExpression(const MatrixOperation<Matrix, Expression2, OPERATIONTYPE_PRODUCT>& op, Matrix& dest)
{
	MatrixOperation<Expression1, Matrix, OPERATIONTYPE_PRODUCT> swapped(op.rightExpr, op.leftExpr);
	expandExpression<state>(swapped, dest);
}


// ### Vector operations ###

//Vector-Vector sum, assign the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_UNINITIALIZED>(const MatrixOperation<Vector, Vector, OPERATIONTYPE_SUM>& op, Vector& dest)
{
	const Vector& v1 = op.leftExpr;
	const Vector& v2 = op.rightExpr;

	mkl_domatadd('R', 'N', 'N', 1, v1.dim, op.leftCoefficient, v1.elements, v1.dim, op.rightCoefficient, v2.elements, v2.dim, dest.elements, v1.dim);
	dest.dim = v1.dim;
}
//Vector-Vector sum, add the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_INITIALIZED>(const MatrixOperation<Vector, Vector, OPERATIONTYPE_SUM>& op, Vector& dest)
{
	dest += op.leftExpr;
	dest += op.rightExpr;
}

//Matrix-Vector product, assign the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_UNINITIALIZED>(const MatrixOperation<Matrix, Vector, OPERATIONTYPE_PRODUCT>& op, Vector& dest)
{
	assert(dest.elements != 0);

	const Matrix& m = op.leftExpr;
	const Vector& v = op.rightExpr;
	Dim mDim = ((op.leftFlag == M_TRANSPOSE) ? Dim{m.columns, m.rows} : Dim{m.rows, m.columns});

	assert(mDim.columns == v.dim);

	cblas_dgemv(CblasRowMajor, (CBLAS_TRANSPOSE)cblasFlagTable[(int)op.leftFlag], m.rows, m.columns, 
					op.leftCoefficient*op.rightCoefficient, m.elements, m.columns, v.elements, 1, 0.0, dest.elements, 1);
	dest.dim = mDim.rows;
}
//Matrix-Vector product, add the result to dest
template<>
inline void expandExpression<EXPRESSIONSTATE_INITIALIZED>(const MatrixOperation<Matrix, Vector, OPERATIONTYPE_PRODUCT>& op, Vector& dest)
{
	assert(dest.elements != 0);

	const Matrix& m = op.leftExpr;
	const Vector& v = op.rightExpr;
	Dim mDim = ((op.leftFlag == M_TRANSPOSE) ? Dim{m.columns, m.rows} : Dim{m.rows, m.columns});

	assert((m.columns == v.dim) && (dest.dim == mDim.rows));

	cblas_dgemv(CblasRowMajor, (CBLAS_TRANSPOSE)cblasFlagTable[(int)op.leftFlag], m.rows, m.columns, 
					op.leftCoefficient*op.rightCoefficient, m.elements, m.columns, v.elements, 1, 1.0, dest.elements, 1);
}

//Scalar-Vector product, assign or add the result to dest
template <ExpressionState state, typename Expression1, typename = std::enable_if_t<std::is_scalar_v<Expression1>>>
inline void expandExpression(const MatrixOperation<Expression1, Vector, OPERATIONTYPE_PRODUCT>& op, Vector& dest)
{
	assert(dest.elements != 0);
	if constexpr(state == EXPRESSIONSTATE_UNINITIALIZED)
	{
		cblas_dscal(dest.dim, op.leftExpr, op.rightExpr.elements, 1);
		dest.dim = op.rightExpr.dim;
	}
	else
	{
		assert(dest.dim == op.rightExpr.dim);
		cblas_daxpy(dest.dim, op.leftExpr, op.rightExpr.elements, 1, dest.elements, 1);
	}
}

template <ExpressionState state, typename Expression2, typename = std::enable_if_t<std::is_scalar_v<Expression2>>>
inline void expandExpression(const MatrixOperation<Vector, Expression2, OPERATIONTYPE_PRODUCT>& op, Vector& dest)
{
	MatrixOperation<Expression2, Vector, OPERATIONTYPE_PRODUCT> swapped(op.rightExpr, op.leftExpr);
	expandExpression<state>(swapped, dest);
}


MatrixConversion::operator Vector()
{
	Vector rowVector(matrix->columns, matrix->elements + index*matrix->columns);
	return rowVector;
}
	
inline Vector& operator*=(r64 a, Vector& v)
{
	cblas_daxpy(v.dim, a - 1.0, v.elements, 1, v.elements, 1);
	return v;
}

inline Vector& operator*=(Vector& v, r64 a)
{
	cblas_daxpy(v.dim, a - 1.0, v.elements, 1, v.elements, 1);
	return v;
}

inline Vector& operator-=(Vector& v1, Vector v2)
{
	assert(v1.dim == v2.dim);

	cblas_daxpy(v1.dim, -1.0, v2.elements, 1, v1.elements, 1);

	return v1;
}

inline Matrix& operator-=(Matrix& m1, Matrix m2)
{
	assert((m1.rows == m2.rows) && (m1.columns == m2.columns));

	cblas_daxpy(m1.rows*m1.columns, -1.0, m2.elements, 1, m1.elements, 1);

	return m1;
}

inline r64 norm(Vector v)
{
	r64 n = cblas_dnrm2(v.dim, v.elements, 1);

	return n;
}

inline void normalize(Vector& v)
{
	r64 n = norm(v);
	v *= 1.0/n;
}

// ####################################################################

//For scalar-matrix multiplication and matrix-matrix addition 
inline void scaleAndAddMatrices(Matrix& result, Matrix m1, Matrix m2, r64 m1Coefficient, r64 m2Coefficient, 
										  MATRIX_OP_FLAG m1Transpose, MATRIX_OP_FLAG m2Transpose)
{
	Dim m1Dim = ((m1Transpose == M_TRANSPOSE) ? Dim{m1.columns, m1.rows} : Dim{m1.rows, m1.columns});
	Dim m2Dim = ((m2Transpose == M_TRANSPOSE) ? Dim{m2.columns, m2.rows} : Dim{m2.rows, m2.columns});

	assert((m1Dim.rows == m2Dim.rows) && (m1Dim.columns == m2Dim.columns));

	result.rows = m1Dim.rows;
	result.columns = m1Dim.columns;

	//NOTE: When transposition etc. is not needed, could also use daxpy
	mkl_domatadd('R', mklFlagTable[(int)m1Transpose], mklFlagTable[(int)m2Transpose], m1Dim.rows, m1Dim.columns,
					 m1Coefficient, m1.elements, m1Dim.columns, m2Coefficient, m2.elements, m2Dim.columns, 
					 result.elements, result.columns);
}

inline void scaleAndAddMatrices(Matrix& result, Matrix m, r64 mCoefficient, MATRIX_OP_FLAG mTranspose)
{
	Dim mDim = ((mTranspose == M_TRANSPOSE) ? Dim{m.columns, m.rows} : Dim{m.rows, m.columns});

	result.rows = mDim.rows;
	result.columns = mDim.columns;

	//NOTE: When transposition etc. is not needed, could also use daxpy
	mkl_domatadd('R', mklFlagTable[(int)mTranspose], 'N', mDim.rows, mDim.columns,
					 mCoefficient, m.elements, mDim.columns, 0.0, m.elements, mDim.columns, 
					 result.elements, result.columns);
}

//For (scalar-)matrix-matrix multiplication (adds the result to 'result')
inline void scaleAndMultiplyMatrices(Matrix& result, Matrix m1, Matrix m2, r64 m1m2Coefficient, r64 resultCoefficient, 
                				  			 MATRIX_OP_FLAG m1Transpose, MATRIX_OP_FLAG m2Transpose)
{
	Dim m1Dim = ((m1Transpose == M_TRANSPOSE) ? Dim{m1.columns, m1.rows} : Dim{m1.rows, m1.columns});
	Dim m2Dim = ((m2Transpose == M_TRANSPOSE) ? Dim{m2.columns, m2.rows} : Dim{m2.rows, m2.columns});

	assert((m1Dim.columns == m2Dim.rows) && (m1Dim.rows == result.rows) && (m2Dim.columns == result.columns));
	//result.rows = m1Dim.rows;
	//result.columns = m2Dim.columns;

	cblas_dgemm(CblasRowMajor, (CBLAS_TRANSPOSE)cblasFlagTable[(int)m1Transpose], (CBLAS_TRANSPOSE)cblasFlagTable[(int)m2Transpose], 
					m1Dim.rows, m2Dim.columns, m1Dim.columns, m1m2Coefficient, m1.elements, m1.columns, 
					m2.elements, m2.columns, resultCoefficient, result.elements, result.columns);
}

inline void scaleAndMultiplyMatrices(Matrix& result, Matrix m1, Matrix m2, r64 m1m2Coefficient, 
												 MATRIX_OP_FLAG m1Transpose, MATRIX_OP_FLAG m2Transpose)
{
	Dim m1Dim = ((m1Transpose == M_TRANSPOSE) ? Dim{m1.columns, m1.rows} : Dim{m1.rows, m1.columns});
	Dim m2Dim = ((m2Transpose == M_TRANSPOSE) ? Dim{m2.columns, m2.rows} : Dim{m2.rows, m2.columns});

	assert((m1Dim.columns == m2Dim.rows)); //&& (m1Dim.rows == result.rows) && (m2Dim.columns == result.columns));
	result.rows = m1Dim.rows;
	result.columns = m2Dim.columns;

	cblas_dgemm(CblasRowMajor, (CBLAS_TRANSPOSE)cblasFlagTable[(int)m1Transpose], (CBLAS_TRANSPOSE)cblasFlagTable[(int)m2Transpose], 
					m1Dim.rows, m2Dim.columns, m1Dim.columns, m1m2Coefficient, m1.elements, m1.columns, //m1Dim.columns and m2Dim.columns ?
					m2.elements, m2.columns, 0.0, result.elements, result.columns);
}


//One iteration of Gram-Schmidt
void orthonormalizeVectorAgainst(Vector& v, Matrix basis)
{
	for(int i = 0; i < basis.rows; i++)
	{
		Vector basisVector(basis(i));
		r64 dot = v*basisVector;
		r64 iNorm = 1.0/norm(basisVector);

		v += -(1.0/norm(basisVector))*dot*basisVector;
	}
 	
	normalize(v);
}

// ### Debug utilities ###

#define str(s) _str(s)
#define _str(s) #s

#ifdef CI_DEBUG
	#define printMatrix(m, rowMajor) _printMatrix(str(m), m, rowMajor)
	#define printPacked(m, rowMajor) _printPacked(str(m), m, rowMajor)
#else
	#define printMatrix(m, rowMajor)
	#define printPacked(m, rowMajor)
#endif


void _printMatrix(char* name, Matrix m, bool rowMajor)
{
	printf("\n\n%s\n", name);
	if(rowMajor)
	{
		for(int i = 0; i < 9*m.columns; i++) {printf("-");}
		printf("\n");
		for(int row = 0; row < m.rows; row++)
		{
			printf("|");
			for(int column = 0; column < m.columns; column++)
			{
				printf(" %8.*g", 4, m(row, column));
			}
			printf("|\n");
		}
		for(int i = 0; i < 9*m.columns; i++) {printf("-");}
		printf("|\n");
	}
	else
	{
		for(int i = 0; i < 9*m.rows; i++) {printf("-\n");}
		printf("\n");
		for(int column = 0; column < m.columns; column++)
		{
			printf("|");
			for(int row = 0; row < m.rows; row++)
			{
				printf(" %8.*g", 4, m(row, column));
			}
			printf(" |\n");
		}
		for(int i = 0; i < 9*m.rows; i++) {printf("-");}
		printf("|\n");
	}

}

void _printMatrix(char* name, Vector v, bool rowMajor)
{
	Matrix matricized(1, v.dim, v.elements);
	_printMatrix(name, matricized, rowMajor);
}

void _printPacked(char* name, Matrix m, bool rowMajor) //Always lower triangular and square
{
	assert(m.rows == m.columns);

	printf("\n\n%s\n", name);
	if(rowMajor)
	{
		for(int i = 0; i < 9*m.columns; i++) {printf("-");}
		printf("\n");
		size_t flatIndex = 0;
		for(int row = 0; row < m.rows; row++)
		{
			printf("|");
			for(int column = 0; column < row+1; column++)
			{
				printf(" %8.*g", 4, m[flatIndex++]);
			}

			for(int column = row+1; column < m.columns; column++)
			{
				printf(" %8.*s", 4, "*");
			}
			printf("|\n");
		}
		for(int i = 0; i < 9*m.columns; i++) {printf("-");}
		printf("|\n");
	}
	else
	{
	}

}

void sort(int* array, size_t length)
{
	for(int j = 1; j < length; j++)
	{
		int compare = array[j];
		int i = j - 1;
		while(i >= 0 && array[i] > compare)
		{
			array[i + 1] = array[i];
			i--;
		}
		array[i + 1] = compare;
	}
}



#endif
