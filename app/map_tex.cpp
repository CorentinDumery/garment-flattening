#include <igl/readOBJ.h>
#include <igl/opengl/glfw/map_texture.cpp>
//#include "tutorial_shared_path.h"
#include "../../external/stb_image/igl_stb_image.cpp"

Eigen::MatrixXd U;
Eigen::MatrixXd V;
Eigen::MatrixXi F;

#define TUTORIAL_SHARED_PATH "../shared"

int main(int argc, char *argv[])
{
  // Load V,U and F from OBJ. F is the same in both files.
  igl::readOBJ(TUTORIAL_SHARED_PATH "/atlas.obj", V, F);
  igl::readOBJ(TUTORIAL_SHARED_PATH "/deformed_atlas.obj", U, F);

  int width, height, bpp, nc =3;
  std::vector<unsigned char> out_rgb;
  //Read image with stbi_load
  unsigned char * in_rgb  = igl::stbi_load(TUTORIAL_SHARED_PATH "/streched_atlas.png", &width, &height, &bpp, nc );

  // in_rgb is now three bytes per pixel, width*height size. Or NULL if load failed.
  if(in_rgb == NULL){
	  printf("image reading failed!\n");
	  return 1;
  }

  //Calling map texture
  printf("Calling map_texture\n");
  igl::opengl::glfw::map_texture(V,F,U,in_rgb,width,height,nc,out_rgb);
  printf("Finished map_texture\n");
  //Write output image with stbi_write
  igl::stbi_write_png(TUTORIAL_SHARED_PATH "/deformed_atlas.png", width, height, bpp, out_rgb.data(), 0);

  igl::stbi_image_free(in_rgb);
}