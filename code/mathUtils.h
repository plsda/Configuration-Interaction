#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <mkl.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef int32_t s32;
typedef double r64;
typedef uint32_t b32;

#define PI64 3.14159265358979323846 
#define IPI64 0.318309886183790671538 
#define IS_EVEN(n) (((n) % 2) == 0)
#define MIN(a,b) ((a) > (b) ? (b) : (a))

#define assert(expr) if(!(expr)){*(int*)0 = 0;}
#define count(array) (sizeof(array)/sizeof((array)[0]))

struct Matrix;
struct Vector;

inline Vector vmmul(Vector v, Matrix m);
inline void mmul(Matrix& result, Matrix m1, Matrix m2, bool, bool);
inline Matrix mmul(Matrix m1, Matrix m2, bool, bool);
inline void mmul(Matrix& result, Matrix m1, Matrix m2);
inline Matrix mmul(Matrix m1, Matrix m2);
inline int orthonormalizeVectors(Matrix basis, r64* workArea);


struct MatrixVectorProduct 
{
	Matrix* m;
	Vector* v;
};

struct Dim
{
	size_t rows;			
	size_t columns;
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

	Matrix(size_t maxRows, size_t maxColumns) : rows(maxRows), columns(maxColumns)
	{
		elements = (r64*)calloc(maxRows*maxColumns, sizeof(r64));
	}

	Matrix(size_t maxRows, size_t maxColumns, size_t activeRows, size_t activeColumns) : 
		rows(activeRows), columns(activeColumns)
	{
		elements = (r64*)calloc(maxRows*maxColumns, sizeof(r64));
	}

	inline r64* operator[](size_t row)
	{
		return (elements + row*columns);
	}

	inline MatrixVectorProduct operator*(Vector& v)
	{
		MatrixVectorProduct result;
		result.m = this;
		result.v = &v;

		return result;
	}

};

struct Vector
{
	r64* elements;
	size_t dim;

	Vector(){}
	Vector(size_t size) : dim(size)
	{
		elements = (r64*)calloc(size, sizeof(r64));
	}
	Vector(size_t size, r64* memory) : dim(size), elements(memory){}

	inline r64& operator[](size_t index)
	{
		return elements[index];
	}

	inline r64 operator*(Vector v)
	{
		return cblas_ddot(dim, elements, 1, v.elements, 1);
	}

	inline Vector& operator=(Vector v)
	{
		cblas_dcopy(dim, v.elements, 1, elements, 1);
		dim = v.dim;

		return *this;
	}

	inline Vector& operator=(const MatrixVectorProduct& mvp)
	{
		Matrix* m = mvp.m;
		Vector* v = mvp.v;

		assert(m->columns == v->dim);

		cblas_dgemv(CblasRowMajor, CblasNoTrans, m->rows, m->columns, 1.0, 
						m->elements, m->columns, v->elements, 1, 0.0, elements, 1);
		dim = v->dim;

		return *this;
	}

};

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


inline void mmul(Matrix& result, Matrix m1, Matrix m2, bool m1Transpose, bool m2Transpose)
{

	int m1Trans = (int)m1Transpose;
	int m2Trans = (int)m2Transpose;
	assert((m1Trans*m1.rows + (1-m1Trans)*m1.columns) == ((1-m2Trans)*m2.rows + m2Trans*m2.columns));

	result.rows = m1Transpose ? m1.columns : m1.rows;
	result.columns = m2Transpose ? m2.rows : m2.columns;
	cblas_dgemm(CblasRowMajor, m1Transpose ? CblasTrans : CblasNoTrans, m2Transpose ? CblasTrans : CblasNoTrans,
					result.rows, result.columns, m1Transpose ? m1.rows : m1.columns, 1.0, m1.elements,
					m1.columns, m2.elements, m2.columns, 0, result.elements, result.columns);
}

inline Matrix mmul(Matrix m1, Matrix m2, bool m1Transpose, bool m2Transpose)
{
	Matrix result;
	result.rows = m1Transpose ? m1.columns : m1.rows;
	result.columns = m2Transpose ? m2.rows : m2.columns;
	result.elements = (r64*)calloc(result.rows*result.columns, sizeof(r64));

	if(result.elements != 0)
	{
		cblas_dgemm(CblasRowMajor, m1Transpose ? CblasTrans : CblasNoTrans, m2Transpose ? CblasTrans : CblasNoTrans,
						result.rows , result.columns, m1Transpose ? m1.rows : m1.columns, 1.0, m1.elements,
						m1.columns, m2.elements, m2.columns, 0, result.elements, result.columns);
	}
	else
	{
	}

	return result;
}

inline void mmul(Matrix& result, Matrix m1, Matrix m2)
{
	mmul(result, m1, m2, false, false);
}

inline Matrix mmul(Matrix m1, Matrix m2)
{
	return mmul(m1, m2, false, false);
}


//One iteration of Gram-Schmidt
void orthonormalizeVectorAgainst(Vector& v, Matrix basis)
{
	Vector bVector(basis.columns, basis.elements);
	Vector temp1(basis.columns);
	Vector temp2(basis.columns);

	for(int i = 0; i < basis.rows; i++)
	{
		bVector.elements = basis[i];
		temp1 = bVector;
		temp2 = bVector;

		r64 bScale = 1.0/norm(bVector);
		bScale *= temp1*v;
		temp2 *= bScale;
		v -= temp2;
	}
 	
	normalize(v);
}

inline int orthonormalizeVectors(Matrix basis, r64* workArea)
{
	int info = (int)LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'O', 'N', basis.columns, basis.rows, basis.elements, 
											 basis.columns, workArea, 0, basis.columns, 0, basis.rows, workArea);

	return info;
}





#endif
