#include "data.h"

Matrix meanFlattenedFace(Matrix matrix);

class FaceSpace {
private:
    Matrix faceDatabase;
    Matrix meanFace;
    Matrix eigenFaceDatabase;
    long numFaces;
    unsigned int numEigenFaces;
    int numOfTopEigenFaces;

    // For adding and removing faces (will come later)
    unsigned int numFacesAdded;
    unsigned int numFacesRemoved;

public:
    // Constructors
    FaceSpace(Matrix faceMatrix, double percentVariance); // The percentVariance parameter allows the user to select what percent of the variance of faces the eigenfaces should cover
    FaceSpace(std::string& path, double percentVariance);





};
