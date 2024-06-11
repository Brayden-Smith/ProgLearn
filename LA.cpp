#include "LA.h"

void thresholdStabilize(Matrix* matrixToStabilize) {
    for (int i = 0; i < matrixToStabilize->getNumRows(); i++) {
        for (int j = 0; j < matrixToStabilize->getNumCols(); j++) {
            if (abs((*matrixToStabilize)(i,j).getRealPart()) < 1e-9) {
                (*matrixToStabilize)(i,j) = ComplexNum(0, (*matrixToStabilize)(i,j).getImagPart());
            }
            if (abs((*matrixToStabilize)(i,j).getImagPart()) < 1e-9) {
                (*matrixToStabilize)(i,j) = ComplexNum((*matrixToStabilize)(i,j).getRealPart(), 0);
            }

        }
    }
}


Matrix matMul(Matrix* lhs, Matrix* rhs) {
    if (lhs->getNumCols() != rhs->getNumRows()) {
        throw std::invalid_argument("matMul: Invalid matrix dimensions");
    }

    Matrix product(lhs->getNumRows(), rhs->getNumCols());
    for (int i = 0; i < lhs->getNumRows(); i++) {
        for (int j = 0; j < rhs->getNumCols(); j++) {
            ComplexNum sum(0, 0);
            for (int k = 0; k < lhs->getNumCols(); k++) {
                sum = sum + ((*lhs)(i,k) * (*rhs)(k,j));
            }
            product(i,j) = sum;
        }
    }
    thresholdStabilize(&product);
    return product;
}

double realTraceOfMatrix(Matrix* matrixToTrace) { // Only defined for real matrices at the moment
    if (matrixToTrace->getNumRows() != matrixToTrace->getNumCols()) {
        throw std::invalid_argument("realTraceOfMatrix: Trace is only defined for square matrices");
    }

    double runningTrace = 0;
    for (int i = 0; i < matrixToTrace->getNumRows(); i++) {
        double numberToAddToRunningTrace = complexNumToDouble((*matrixToTrace)(i,i));
        runningTrace += numberToAddToRunningTrace;
    }

    return runningTrace;
}

Matrix conjTranspose(Matrix* matrixToTranspose) {
    Matrix matrixToReturn(matrixToTranspose->getNumCols(), matrixToTranspose->getNumRows());
    for (int i = 0; i < matrixToTranspose->getNumRows(); i++) {
        for (int j = 0; j < matrixToTranspose->getNumCols(); j++) {
            ComplexNum transposedEntry = complexConjugate((*matrixToTranspose)(i,j));
            matrixToReturn(j,i) = transposedEntry;
        }
    }
    return matrixToReturn;
}

Matrix transpose(Matrix* matrixToTranspose) {
    Matrix matrixToReturn(matrixToTranspose->getNumCols(), matrixToTranspose->getNumRows());
    for (int i = 0; i < matrixToTranspose->getNumRows(); i++) {
        for (int j = 0; j < matrixToTranspose->getNumCols(); j++) {
            ComplexNum transposedEntry = (*matrixToTranspose)(i,j);
            matrixToReturn(j,i) = transposedEntry;
        }
    }
    return matrixToReturn;
}

double frobeniusNorm(Matrix* matrixToNorm) {
    Matrix conjugateTransposeMatrix = conjTranspose(matrixToNorm);
    Matrix matrixStarMatrix = matMul(&conjugateTransposeMatrix, matrixToNorm);
    //std::cout << "Matrix mul has dimensions " << matrixStarMatrix.numRows << "x" << matrixStarMatrix.numCols << std::endl;
    double trace = realTraceOfMatrix(&matrixStarMatrix);
    return sqrt(trace);
}

ComplexNum InnerProduct(Matrix& u, Matrix& v) {
    if (u.getNumCols() != 1 || v.getNumCols() != 1) {
        throw std::invalid_argument("InnerProduct: Vectors must be nx1 matrices");
    }

    ComplexNum result(0, 0);
    for (int i = 0; i < u.getNumRows(); i++) {
        ComplexNum vEntryConjugate = complexConjugate(v(i, 0));
        result += (u(i, 0) * vEntryConjugate);
    }
    return result;
}


// Method to normalize vectors in a matrix (alters matrix upon which it is called)
void normalizeVectorsInMatrix(Matrix* pointerToMatrix) {
    for (int i = 0; i < pointerToMatrix->getNumCols(); i++) {

        double Norm = 0;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            double magnitude = magnitudeOfNumber((*pointerToMatrix)(j,i));
            Norm += (magnitude * magnitude);
        }
        if (Norm == 1) {
            return;
        }
        Norm = sqrt(Norm);
        double oneOverMagnitude = 1/Norm;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            (*pointerToMatrix)(j,i) = (*pointerToMatrix)(j,i) * oneOverMagnitude;
        }

    }
}

Matrix GramSchmidt(Matrix const& M) {
    Matrix result(M.getNumRows(),M.getNumCols());

    for(int i = 0; i < M.getNumCols(); i++) {
        Matrix Vk = M[i];
        Matrix Uk = Vk;
        for(int j = 0; j < i; j++) {

            Matrix Un = result[j];
            ComplexNum VdotU = InnerProduct(Vk,Un);
            ComplexNum UdotU = InnerProduct(Un,Un);

            //compute the projection
            ComplexNum quotient = (VdotU / UdotU);
            Un = Un * quotient;
            Uk = Uk - Un;
        }
        result.columnAssign(i,&Uk);
    }

    normalizeVectorsInMatrix(&result);
    return result;
}

void swapRowsInMatrix(Matrix* matrixToSwap, int posRowOne, int posRowTwo) {

    Matrix temp(1, matrixToSwap->getNumCols());
    for (int i = 0; i < matrixToSwap->getNumCols(); i++) {
        temp(0, i) = (*matrixToSwap)(posRowOne, i);
    }

    for (int i = 0; i < matrixToSwap->getNumCols(); i++) {
        (*matrixToSwap)(posRowOne, i) = (*matrixToSwap)(posRowTwo, i);
    }

    for (int i = 0; i < matrixToSwap->getNumCols(); i++) {
        (*matrixToSwap)(posRowTwo, i) = temp(0, i);
    }
}

Matrix createAugmentedMatrix(Matrix* matrixA, Matrix* matrixB) {
    if (matrixA->getNumRows() != matrixB->getNumRows()) {
        throw std::invalid_argument("createAugmentedMatrix: Matrices must have same number of rows");
    }
    Matrix matrixToReturn(matrixA->getNumRows(), matrixA->getNumCols() + matrixB->getNumCols());

    for (int i = 0; i < matrixA->getNumRows(); i++) {
        for (int j = 0; j < matrixA->getNumCols(); j++) {
            matrixToReturn(i, j) = (*matrixA)(i, j);
        }
        for (int j = 0; j < matrixB->getNumCols(); j++) {
            matrixToReturn(i, j + matrixA->getNumCols()) = (*matrixB)(i, j);
        }
    }
    return matrixToReturn;
}


// Row reduction to check for linear independence will be a separate method called "row reduction"
Matrix gaussianElimination(Matrix* matrix, Matrix* vector) {

    Matrix augmentedMatrixPrePointer = createAugmentedMatrix(matrix, vector); // Temporary fix
    Matrix* augmentedMatrix = &augmentedMatrixPrePointer;

    //std::cout << "Augmented matrix to perform row reduction on:\n" << augmentedMatrixPrePointer << std::endl;

    if (augmentedMatrix->getNumRows() + 1 != augmentedMatrix->getNumCols()) {
        throw std::invalid_argument("gaussianElimination: Augmented matrix must be of size n x n+1");
    }
    // Forward elimination
    for (int i = 0; i < augmentedMatrix->getNumCols() - 2; i++) {
        // Get the largest element in column and swap rows for numerical stability

        ComplexNum maxElementInColumn(0,0);
        int maxElementInColumnRowPos = i;
        for (int j = i; j < augmentedMatrix->getNumRows(); j++) {
            if (magnitudeOfNumber((*augmentedMatrix)(j, i)) > magnitudeOfNumber(maxElementInColumn)) {
                maxElementInColumnRowPos = j;
                maxElementInColumn = (*augmentedMatrix)(j, i);
            }
        }
        if (magnitudeOfNumber(maxElementInColumn) < 1e-6) {
            throw std::invalid_argument("gaussianElimination: Matrix is singular!");
        }

        if (maxElementInColumnRowPos != i) {
            swapRowsInMatrix(augmentedMatrix, i, maxElementInColumnRowPos);
        }

        //std::cout << "LARGEST ELEMENT IN COLUMN IS " << maxElementInColumn << std::endl;

        // Use operations to ensure the column is consistent with a matrix in REF

        for (int j = i + 1; j < augmentedMatrix->getNumRows(); j++) {
            // Add multiples of j = i to each row from i + 1 to the end of column to get zeros
            ComplexNum scalarMultipleForRowi = ((*augmentedMatrix)(j, i))/(*augmentedMatrix)(i, i);
            scalarMultipleForRowi = scalarMultipleForRowi * -1;

            //std::cout << augmentedMatrixPrePointer << std::endl;

            //std::cout << "scalar multiple for row " << j << " is " << scalarMultipleForRowi << std::endl;

            for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
                (*augmentedMatrix)(j, k) = (*augmentedMatrix)(j, k) + ((*augmentedMatrix)(i, k) * scalarMultipleForRowi);
                //std::cout << (*augmentedMatrix)(j, k) << std::endl;
            }
            //std::cout << augmentedMatrixPrePointer << std::endl;

        }

    }
    //std::cout << "Upper triangular matrix is\n" << augmentedMatrixPrePointer << "\n" << std::endl;

    // Execute the same procedure now considering the entries above the diagonal (this can be replaced with straight back substitution in the future)
    for (int i = 0; i < augmentedMatrix->getNumCols() - 1; i++) {
        // Get the largest element in column and swap rows for numerical stability
        if (magnitudeOfNumber((*augmentedMatrix)(i, i)) < 1e-6) {
            throw std::invalid_argument("gaussianElimination: Matrix is singular!");
        }
        //std::cout << "here?" << std::endl;

        for (int j = i - 1; j >= 0; j--) {

            // Add multiples of j = i to each row from i + 1 to the end of column to get zeros
            ComplexNum scalarMultipleForRowi = ((*augmentedMatrix)(j, i))/(*augmentedMatrix)(i, i);
            scalarMultipleForRowi = scalarMultipleForRowi * -1;

            //std::cout << "scalar multiple for row " << i + 1 << " is " << scalarMultipleForRowi << std::endl;

            for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
                (*augmentedMatrix)(j, k) = (*augmentedMatrix)(j, k) + ((*augmentedMatrix)(i, k) * scalarMultipleForRowi);
            }

        }

    }

    //std::cout << augmentedMatrixPrePointer << std::endl;
    for (int j = 0; j < augmentedMatrix->getNumRows(); j++) {
        ComplexNum reciprocal = ComplexNum(1, 0) / (*augmentedMatrix)(j,j);
        for (int k = 0; k < augmentedMatrix->getNumCols(); k++) {
            (*augmentedMatrix)(j,k) = (*augmentedMatrix)(j,k) * reciprocal;
        }
    }

    Matrix vectorToReturn(vector->getNumRows(), vector->getNumCols());
    for (int j = 0; j < vector->getNumRows(); j++) {
        vectorToReturn(j, 0) = (*augmentedMatrix)(j, augmentedMatrix->getNumCols() - 1);
    }

    return vectorToReturn;
}

Matrix inverseMatrix(Matrix* matrixToInvert) { // Really, really slow. That being said, we can speed this up in the future if needed, but at the moment I want "easy" code for debugging
    if (matrixToInvert->getNumRows() != matrixToInvert->getNumCols()) {
        throw std::invalid_argument("inverseMatrix: Matrix must be nxn!");
    }

    Matrix outputMatrix(matrixToInvert->getNumRows(), matrixToInvert->getNumCols());
    for (int i = 0; i < matrixToInvert->getNumRows(); i++) {
        Matrix unitVector(matrixToInvert->getNumRows(), 1);
        unitVector(i, 0) = 1;
        Matrix resultVector = gaussianElimination(matrixToInvert, &unitVector);

        for (int j = 0; j < matrixToInvert->getNumRows(); j++) {
            outputMatrix(j, i) = resultVector(j , 0);
        }
        //std::cout << outputMatrix[i] << std::endl;
    }

    return outputMatrix;
}

//gets magnitude of column matrix
double VectorNorm(Matrix* pointerToMatrix) {
    double Norm = 0;
    for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
        double magnitude = magnitudeOfNumber((*pointerToMatrix)(j,0));
        Norm += (magnitude * magnitude);
    }
    if(Norm == 1) {
        return Norm;
    }
    return sqrt(Norm);
}

Matrix unitVector(int k, int dim) {
    if (k > dim || dim < 1  || k < 0) {
        throw std::invalid_argument("unitVector: Invalid input");
    }
    Matrix e(dim,1);
    e(k,0) = 1;
    return e;
}

Matrix identityMatrix(int dim) {
    if(dim < 1) {
        throw std::invalid_argument("identityMatrix: Invalid input");
    }
    Matrix I(dim,dim);

    //replace each column with unit vector
    for(int i = 0; i < dim; i++) {
        Matrix ek = unitVector(i,dim);
        I.columnAssign(i,&ek);
    }

    return I;
}

Matrix frontFillVec(Matrix M, int dim, ComplexNum const& Fill) {
    if(M.getNumCols() > 1 || dim < 1 || dim < M.getNumCols()) {
        throw std::invalid_argument("frontFillVec: Invalid input");
    }


    Matrix result(dim,1);
    for(int i = 0; i < dim; i++) {
        //fill with num until getting to matrix with data
        if(i < dim - M.getNumRows()) {
            result(i,0) = Fill;
        }
        else {
            result(i,0) = M(i - (dim - M.getNumRows()),0);
        }
    }
    return result;
}

Matrix minorMatrix(Matrix M,int iDel, int jDel) {
    Matrix Minor(M.getNumRows() -1,M.getNumCols() -1);
    for(int i = 0; i < M.getNumRows(); i++) {
        if (i < iDel) {
            for(int j = 0; j < M.getNumCols(); j++) {{
                if(j < jDel) {
                    Minor(i,j) = M(i,j);
                }
                else if(j > jDel) {
                    Minor(i,j-1) = M(i,j);
                }
            }}
        }
        else if(i > iDel) {
            for(int j = 0; j < M.getNumCols(); j++) {{
                if(j < jDel) {
                    Minor(i-1,j) = M(i,j);
                }
                else if(j > jDel) {
                    Minor(i-1,j-1) = M(i,j);
                }
            }}
        }
    }
    return Minor;
}

Matrix householderReflection(Matrix* x) {
    if (x->getNumCols() > 1) {
        throw std::invalid_argument("householderTransform: Invalid input");
    }

    Matrix y = (*x)(0,0).sign() * VectorNorm(x) * unitVector(0,x->getNumRows());
    y = *x + y;
    return y;
}

std::vector<Matrix> QRDecomp(Matrix const& A) {

    Matrix R = A;
    std::vector<Matrix> HTransforms;
    for(int i = 0; i < A.getNumCols(); i++) {
        if(i > 0) {
            R = minorMatrix(R,0,0);
        }


        //get transformation on column
        Matrix a = R[0];
        Matrix v = householderReflection(&a);
        Matrix vT = conjTranspose(&v);


        //get the transformation on the matrix
        Matrix H = identityMatrix(R.getNumRows());
        H = H - ((ComplexNum(2,0) / InnerProduct(v,v) * (matMul(&v,&vT)) ));
        HTransforms.push_back(H);

        R = matMul(&H,&R);


    }

    for(int i = 1; i < HTransforms.size(); i++) {
        Matrix resized = identityMatrix(HTransforms[0].getNumCols());
        for(int j = HTransforms[0].getNumCols()-HTransforms[i].getNumCols(); j < HTransforms[0].getNumCols(); j++){
            Matrix toFrontLoad = (HTransforms[i])[j-(HTransforms[0].getNumCols()-HTransforms[i].getNumCols())];
            Matrix processed = frontFillVec(toFrontLoad,HTransforms[0].getNumCols(),ComplexNum(0,0));
            resized.columnAssign(j,&processed);
        }
        HTransforms[i] = resized;
    }

    Matrix Q = HTransforms.back();
    for(int i = HTransforms.size() - 2; i >= 0; i--) {
        Q = matMul(&HTransforms[i],&Q);
    }

    R = A;
    for(int i = 0; i < HTransforms.size(); i++) {
        R = matMul(&HTransforms[i],&R);
    }
    return {Q,R};
}


bool isUpperTriangular(Matrix* mat) {
    for (int i = 0; i < mat->getNumCols(); ++i) {
        for (int j = i + 1; j < mat->getNumRows(); ++j) {
            if (magnitudeOfNumber((*mat)(j, i)) > 1e-9) {
                return false;
            }
        }
    }
    return true;
}

ComplexNum mu(Matrix* matrix) {
    ComplexNum delta = ((*matrix)(0, 0) - (*matrix)(1, 1))/2;
    ComplexNum radicalMult = sqrt((delta * delta) + ((*matrix)(0, 1)) * ((*matrix)(1, 0)));

    double sign;
    if (delta.getRealPart() >= 0) {
        sign = 1;
    } else {
        sign = -1;
    }

    ComplexNum WilkinsonShift = delta - ((sign) * radicalMult);

    return WilkinsonShift;
}

std::vector<ComplexNum> eigenvalues(Matrix* matrix) {
    if (matrix->getNumRows() != matrix->getNumCols()) {
        throw std::invalid_argument("eigenvalues: matrix must be n x n!");
    }

    if (matrix->getNumCols() == 1 && matrix->getNumRows() == 1) { // Trivial case
        std::vector<ComplexNum> ret;
        ret.push_back(ComplexNum((*matrix)(0, 0)));
        return ret;
    }

    Matrix matrixToIterate = (*matrix);
    // Check to see if matrix is already upper triangular


    bool akIsUpperTriangular = isUpperTriangular(&matrixToIterate);


    while (!akIsUpperTriangular) {
        std::cout << "Matrix to iterate:\n" << matrixToIterate << std::endl;
        std::vector<Matrix> currentQR = QRDecomp(matrixToIterate);
        std::cout << "Q is:\n" << currentQR[0] << std::endl;
        std::cout << "R is:\n" << currentQR[1] << std::endl;
        matrixToIterate = matMul(&currentQR[1], &currentQR[0]);
        //std::cout << matrixToIterate << std::endl;

        // Check if matrixToIterate is upper triangular
        akIsUpperTriangular = isUpperTriangular(&matrixToIterate);
    }
    //std::cout << "akIsUpperTriangular and ak is:\n" << matrixToIterate << std::endl;

    std::set<ComplexNum> setOfEigenvalues;
    for (int i = 0; i < matrixToIterate.getNumRows(); i++) {
        double numImagPart = (matrixToIterate(i, i)).getImagPart();
        double numRealPart = (matrixToIterate(i, i)).getRealPart();
        ComplexNum newCNum(numRealPart, numImagPart);
        setOfEigenvalues.insert(newCNum);
    }

    std::vector<ComplexNum> vectorOfEigenvalues;
    for (const ComplexNum& element : setOfEigenvalues) {
        double numImagPart = element.getImagPart();
        double numRealPart = element.getRealPart();
        ComplexNum newCNum(numRealPart, numImagPart);
        vectorOfEigenvalues.push_back(newCNum);
    }

    return vectorOfEigenvalues;

}

//assumes all values in the matrix have uniform weight (add another functions that uses weight vector??)
ComplexNum expectedValue(Matrix* Z) {
    if(Z->getNumCols() > 1) {
        throw std::invalid_argument("expectedValue: Not a complex random variable");
    }

    ComplexNum E(0,0);
    for(int i = 0; i < Z->getNumRows(); i++) {
        E += (*Z)(i,0);
    }
    E = E * (1.0/(Z->getNumRows()));
    return E;
}

ComplexNum covariance(Matrix* Z, Matrix* W) {
    Matrix conjW = W->conjugate();
    Matrix ZconjW = matMul(Z,&conjW);
    ComplexNum EZconjW = expectedValue(&ZconjW);

    ComplexNum EZ = expectedValue(Z);
    ComplexNum EW = expectedValue(W);

    return EZconjW - (EZ * EW);
}

Matrix covarianceMatrix(Matrix* M) {

}
