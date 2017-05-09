#include <stdio.h>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

typedef unsigned char byte;

vector<byte> matToBytes(Mat image)
{
    int size = image.total() * image.elemSize();
    vector<byte> img_bytes(size);
    img_bytes.assign(image.datastart, image.dataend);
    return img_bytes;
}

Mat bytesToMat(vector<byte> bytes,int width,int height)
{
    Mat image = Mat(height,width,CV_8UC3,bytes.data()).clone(); // make a copy
    return image;
}

int main(int argc, char** argv )
{
    if ( argc != 2 )
    {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }

    Mat image;
    image = imread( argv[1], 1 );

    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }
    
    Size sz = image.size(); //used to find image width and height for bytesToMat
    int imageWidth = sz.width;
    int imageHeight = sz.height;
    cout << imageHeight << endl;
    vector<byte> img_bytes = matToBytes(image);
    Mat newimage = bytesToMat(img_bytes, imageWidth, imageHeight);
    imwrite("newimage.jpeg", newimage);
    namedWindow("Display Image", WINDOW_AUTOSIZE );
    imshow("Display Image", newimage);

    return 0;
}
