//! Opens .tga -pictures
#include <iostream>
#include <fstream>
#include <memory.h>

#ifndef __TGA_HPP
#define __TGA_HPP

#define IMG_OK				0x1
#define IMG_ERR_NO_FILE		0x2
#define IMG_ERR_MEM_FAIL	0x4
#define IMG_ERR_BAD_FORMAT	0x8
#define IMG_ERR_UNSUPPORTED 0x40

class TGA {
public:
    //! Constructor
	TGA();
	//! Destructor
	~TGA();
	//! Loads a file of type .tga
	int Load(char* szFilename);
	//! Returns the bits/pixel
	int GetBPP();
	//! Returns the picture width
	int GetWidth();
	//! Returns the picture heigth
	int GetHeight();
	//! Returns a pointer to image data
	unsigned char* GetImg();       // Return a pointer to image data
	//! Returns a pointer to VGA palette
	unsigned char* GetPalette();   // Return a pointer to VGA palette

private:
    //! Width of the image
	short int iWidth;
	//! Height of the image
	short int iHeight;
	//! Bits/pixel
	short int iBPP;
	//! Image Size in bytes
	unsigned long lImageSize;
	//! Encoding of pictures
	char bEnc;
	//! pointer to image data
	unsigned char *pImage;
	//! pointer to palette data
	unsigned char *pPalette;
	//! pointer to raw data
	unsigned char *pData;

	// Internal workers
	//! Read file-header
	int ReadHeader();
	//! Load raw data
	int LoadRawData();
	//! Load the RLE compressed data
	int LoadTgaRLEData();
	//! Load the palette
	int LoadTgaPalette();
	//! Convert BGR to RGB color
	void BGRtoRGB();
	//! Flip the image (vertically)
	void FlipImg();
 };

#endif
