/*
 * =====================================================================================
 *
 *       Filename:  savepng.c
 *
 *    Description:  Grabs pixels from OpenGL buffer and saves them to RGBA32 png
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California Davis
 *
 * =====================================================================================
 */

#include "savepng.h"
#include <png.h>
#include <stdio.h>
#include <GL/gl.h>
#include <stdlib.h>

/*
 * Saves the current OpenGL buffer to a .png in the ./pngs subdirectory of the
 * current working directory.  As an example If filename points to the char
 * array "test\0", and k == 5, the filename written will be of the form:
 *
 * ./pngs/test0005.png
 *
 */
int SavePNG(int width, int height, char *basefilename, int k)
{
  int i;
  png_bytep *row_pointers;
  unsigned char *buffer;
  char filename[50];
  png_structp png_ptr;
  png_infop info_ptr;
  FILE *fp;
  sprintf(filename, "%s%04d.png", basefilename, k);

  // Open file for writing (binary mode)
  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "capture: Couln't open output file \"%s\"\n", filename);
    return 1;
  }

  /*  Initialize PNG structs */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); 
  if (!png_ptr) {
    fprintf(stderr, "capture: Can't initialize png_ptr");
    return 1;
  }
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
     png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
     fprintf(stderr, "capture: Can't initialze info_ptr");
     return 1;
  }

  /*  Initialize PNG error jump */
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    fprintf(stderr, "capture: Unknown error");
    return 1;
  }

  /*  Give PNG the file handle */
  png_init_io(png_ptr, fp);

  /*  Set PNG options/info.  Assumes 24 bit color buffer, change this
   *  if you need to.
   */
  png_set_IHDR(png_ptr, info_ptr, 
      width,                          /* Width */
      height,                         /* Height */
      8,                              /* Bit depth */ 
      PNG_COLOR_TYPE_RGB_ALPHA,       /* Color type */
      PNG_INTERLACE_NONE,             /* Interlacing */
      PNG_COMPRESSION_TYPE_DEFAULT,   /* Compression */
      PNG_FILTER_TYPE_DEFAULT);       /* Filter method */

  /*  Set up row pointers.  OpenGL stores the buffer in reverse
   *  row order to PNG (bottom-to-top instead of top-to-bottom),
   *  so we flip them here.
   */
  buffer = (unsigned char *) malloc(width * height * 4);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0,0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer);

  row_pointers = png_malloc(png_ptr, height * sizeof(png_bytep));
  for (i = 0; i< height; i++)
    row_pointers[i] = &buffer[(height - i - 1) * width* 4];
  png_set_rows(png_ptr, info_ptr, row_pointers);

  /*  Write the PNG */
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

  /*  Free up */
  png_destroy_write_struct(&png_ptr, &info_ptr);
  free(buffer);
  fclose(fp);
  return 0;
} // 
