#include "data.h"

Matrix meanFlattenedFace(Matrix matrix);

class FaceSpace {
private:
    Matrix faceDatabase;
    Matrix meanFace;
    Matrix eigenFaceDatabase;
    long numFaces;
    unsigned int numEigenFaces;
    unsigned int numOfTopEigenFaces;

    // For adding and removing faces (will come later)
    unsigned int numFacesAdded;
    unsigned int numFacesRemoved;

public:
    // Constructors
    FaceSpace(Matrix faceMatrix, double percentVariance); // The percentVariance parameter allows the user to select what percent of the variance of faces the eigenfaces should cover
    FaceSpace(std::string& path, double percentVariance);

    long getNumFaces();
    unsigned int getNumEigenFaces();
    unsigned int getNumTopEigenFaces();

    unsigned int getNumFacesAdded();
    unsigned int getNumFacesRemoved();

    // todo dislay face, match face, remove face, and add face methods

};
