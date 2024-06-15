#pragma once
#include "Matrix.h"
#include <fstream>

class pngImage {
private:
    Matrix pictureData;
    bool isGrayScale;

    //height and width are 4 bytes so int is nice
    int width;
    int height;

    //1 byte
    char bitDepth;
    char colorType;
    char compressionMethod;
    char filterMethod;
    char interlaceMethod;


public:
    pngImage(std::string filePath);
    void convertGray();
};
