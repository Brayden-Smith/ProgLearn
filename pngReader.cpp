#include "pngReader.h"

Matrix png2MonoMatrix(const std::string& filePath) {
    int width, height, channels;
    //file path has to be a c style string
    unsigned char *img = stbi_load(filePath.c_str(),&width,&height,&channels,0);

    if (!img) {
        throw std::invalid_argument("pngReader: Could not load image");
    }

    Matrix result(height,width);

    for (int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            //library flattens the image so we need to find what we want
            int pixelIndex = (i * width + j) * channels;

            unsigned char r = img[pixelIndex];
            unsigned char g = img[pixelIndex + 1];
            unsigned char b = img[pixelIndex + 2];
            unsigned char grayscale = static_cast<unsigned char>((r + g + b) / 3);

            result(i,j) = int(grayscale);
        }
    }

    stbi_image_free(img);
    return result;
}