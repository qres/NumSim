/*
Copyright (C) 2015   Malte Brunn

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//------------------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tga.hpp"

int bye() {
  printf("\nSo long, and thanks for all the fish...\n");
  return 0;
}

int getOpt(const char *opt, int &argc, char **&argv) {
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], opt) == 0) {
      return i;
    }
  }
  return 0;
}

typedef struct {
  double re;
  double omg;
  double alpha;
  double dt;
  double tend;
  double eps;
  double tau;
  int iter;
} param_t;

typedef struct {
  int s_x;
  int s_y;
  double l_x;
  double l_y;
  double v_x;
  double v_y;
  double p;
  int pre;
} geom_t;

int main(int argc, char **argv) {
  int pos;
  int cnt = 0;
  int cnt2 = 0;
  char filename[100] = "earth";
  char *c_map;
  c_map = NULL;
  param_t param;
  param.alpha = 0.9;
  param.re = 1000;
  param.omg = 1.7;
  param.dt = 0;
  param.tend = 50;
  param.eps = 0.001;
  param.tau = 0.5;
  param.iter = 100;
  geom_t geom;
  geom.l_x = 1;
  geom.l_y = 1;
  geom.p = 0;
  geom.s_x = 128;
  geom.s_y = 128;
  geom.v_x = 1;
  geom.v_y = 0;
  geom.pre = 0;

  pos = getOpt("-h", argc, argv);
  if (!pos)
    pos = getOpt("--help", argc, argv);
  if (pos) {
    printf("Magrathea - Custom-Made, Luxury Planet Building.\n\n");
    printf("Usage: magrathea [arguments]\tGenerates a custom-made, luxury "
           "world\n\n");
    printf("General arguments:\nArguments must be seperated by a space from "
           "their additional parameters.\n");
    printf("    -h, --help\t\tPrints this message.\n");
    printf("    -o <name>\t\tSaves the world and its parameters in "
           "\"<name>.param\"\n\t\t\tand \"<name>.geom\". The default filenames "
           "are\n\t\t\t\"earth.param\" and \"earth.geom\".\n");
    printf("    -v, --version\tPrints the version. Or the meaning of live and "
           "stuff...\n\t\t\tTry it yourself.\n");
    printf("\n\nParameters:\nIf some of these arguments are given, the "
           "\".param\" file is written.\n");
    printf("    -alpha <float>\tSets the mysterious alpha value to <float>. "
           "Default: 0.9\n");
    printf("    -dt <float>\t\tSomething with time. It uses <float> as value "
           "though.\n\t\t\tDefault: 0 ...Almost forgot to tell you 0 is "
           "somehow,\n\t\t\tum well, special.\n");
    printf("    -eps <float>\tTolerance or so... Usually very very low. "
           "Default: 0.001\n");
    printf("    -iter <int>\t\tSlows down the planet, but I was told this is "
           "good.\n\t\t\tDefault: 100\n");
    printf("    -omg <float>\tAnother spooky parameter. Its default is 1.7 "
           "though.\n");
    printf("    -re <float>\t\tDefines fluid dynamics. Default: 1000\n");
    printf("    -tau <float>\tIt does things... Planets are complicated "
           "dude.\n\t\t\tIt's set to 0.5 and I wouldn't screw with it.\n");
    printf("    -tend <float>\tSelf-destruction time. Default: 50.0\n");
    printf("\n\nGeometry:\nIf some of these arguments are given, the \".geom\" "
           "file is written.\nThe default world is a driven cavity. "
           "Pre-defined worlds must be set before\nchanging other "
           "parameters.\n");
    printf("    -pre <num>\t\tChoose among some ready-made worlds.\n\t\t\t  1: "
           "Simple channel\n\t\t\t  2: Pressure driven channel\n\t\t\t  3: "
           "Channel with step\n\t\t\t  4: Channel with barrier\n\t\t\tThese "
           "worlds might change the default values.\n");
    printf("    -pd1 <float>\tTweak your ready-made worlds for your needs.\n\t\t\t"
           "Be aware, that this parameter might be ignored for\n\t\t\t"
           "some worlds. Some say it behaves like 'b/2'.\n");
    printf("    -pd2 <float>\tInsane!!! Even more personalization. Although it is\n"
           "\t\t\ta confusing one. This is called alpha as well. (Saay what?!?!).\n");
    printf("    -load <file>\tLoad a TGA picture as geometry. Sounds "
           "impossible,\n\t\t\tbut it works.\n");
    printf("    -length <x>x<y>\tReal real-world world size. Let it set to "
           "1.0x1.0\n\t\t\tunless you do fancy stuff.\n");
    printf("    -pressure <p>\t...this is our last dance, this is "
           "ourselves\n\t\t\tunder pressure, under pressure... Oh, "
           "nevermind.\n\t\t\tIt is 0.0.\n");
    printf("    -size <x>x<y>\tSize matters! Default is 128x128 and that "
           "really\n\t\t\tis a huge world man.\n");
    printf("    -speed <x>x<y>\tInital speed of light... Just kidding, only "
           "speed of\n\t\t\tfluids. Some people care, but I think its fine "
           "at\n\t\t\t1.0x0.0. Whatever that means.\n");
    return bye();
  }
  pos = getOpt("-v", argc, argv);
  if (!pos) {
    pos = getOpt("--version", argc, argv);
  }
  if (pos) {
    printf("42\n");
    return 0;
  }
  pos = getOpt("-o", argc, argv);
  if (pos) {
    sscanf(argv[pos + 1], "%s", filename);
  }

  pos = getOpt("-alpha", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.alpha);
  }

  pos = getOpt("-dt", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.dt);
  }

  pos = getOpt("-eps", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.eps);
  }

  pos = getOpt("-iter", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%i", &param.iter);
  }

  pos = getOpt("-omg", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.omg);
  }

  pos = getOpt("-re", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.re);
  }

  pos = getOpt("-tau", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.tau);
  }

  pos = getOpt("-tend", argc, argv);
  if (pos) {
    cnt++;
    sscanf(argv[pos + 1], "%lf", &param.tend);
  }

  if (cnt) {
    FILE *handle;
    char name[100];
    sprintf(name, "%s.param", filename);
    handle = fopen(name, "w");
    fprintf(handle, "re = %lf\n", param.re);
    fprintf(handle, "omg = %lf\n", param.omg);
    fprintf(handle, "alpha = %lf\n", param.alpha);
    fprintf(handle, "dt = %lf\n", param.dt);
    fprintf(handle, "tend = %lf\n", param.tend);
    fprintf(handle, "iter = %i\n", param.iter);
    fprintf(handle, "eps = %lf\n", param.eps);
    fprintf(handle, "tau = %lf\n", param.tau);
    fclose(handle);
  }

  pos = getOpt("-pre", argc, argv);
  if (pos) {
    cnt2++;
    sscanf(argv[pos + 1], "%i", &geom.pre);
  }
  
  int pre_data1 = 0;
  float pre_b = 0;
  pos = getOpt("-pd1", argc, argv);
  if (pos) {
    sscanf(argv[pos + 1], "%f", &pre_b);
  }

  float pre_a = 45;
  pos = getOpt("-pd2", argc, argv);
  if (pos) {
    sscanf(argv[pos + 1], "%f", &pre_a);
  }
  if (pre_a < 0 || pre_a > 90) pre_a = 45;
  pre_a = pre_a/180.0*M_PI;

  switch (geom.pre) {
  case 1:
    geom.l_x = 5;
    geom.s_x = 160;
    geom.s_y = 32;
    break;
  case 2:
  case 3:
  case 4:
    geom.l_x = 5;
    geom.s_x = 160;
    geom.s_y = 32;
    geom.v_x = 0;
    geom.p = 0.1;
    break;
  default:
    break;
  };

  pos = getOpt("-length", argc, argv);
  if (pos) {
    cnt2++;
    sscanf(argv[pos + 1], "%lfx%lf", &geom.l_x, &geom.l_y);
  }

  pos = getOpt("-size", argc, argv);
  if (pos) {
    cnt2++;
    sscanf(argv[pos + 1], "%ix%i", &geom.s_x, &geom.s_y);
  }

  pos = getOpt("-speed", argc, argv);
  if (pos) {
    cnt2++;
    sscanf(argv[pos + 1], "%lfx%lf", &geom.v_x, &geom.v_y);
  }

  pos = getOpt("-pressure", argc, argv);
  if (pos) {
    cnt2++;
    sscanf(argv[pos + 1], "%lf", &geom.p);
  }
  
  if (pre_b == 0)
    pre_b = geom.l_y * 0.5;
  pre_data1 = geom.s_y * pre_b/geom.l_y;

  pos = getOpt("-load", argc, argv);
  if (pos) {
    cnt2++;
    geom.pre = 5;
    TGA tga;
    tga.Load(argv[pos + 1]);
    if (tga.GetBPP() != 24) {
      printf("Only 24bpp TGA images are supported.\n");
      return bye();
    }
    geom.s_x = tga.GetWidth();
    geom.s_y = tga.GetHeight();
    geom.l_y = geom.l_x * double(geom.s_y) / geom.s_x;
    c_map = new char[geom.s_x * geom.s_y];
    for (int i = 0; i < geom.s_x * geom.s_y; ++i) {
      unsigned int data = tga.GetImg()[3 * i] |
                          (tga.GetImg()[(3 * i) + 1] << 8) |
                          (tga.GetImg()[(3 * i) + 2] << 16);
      switch (data) {
      case 0x00000000:
        c_map[i] = '#';
        break;
      case 0x000000FF:
        c_map[i] = 'V';
        break;
      case 0x0000FFFF:
        c_map[i] = 'H';
        break;
      case 0x0000FF00:
        c_map[i] = '|';
        break;
      case 0x00FFFF00:
        c_map[i] = '-';
        break;
      case 0x00FF0000:
        c_map[i] = 'O';
        break;
      case 0x00FF00FF:
        c_map[i] = 'I';
        break;
      default:
        c_map[i] = ' ';
        break;
      };
    }
    for (int y = 0; y < geom.s_y; ++y) {
      for (int x = 0; x < geom.s_x; ++x) {
        int check = 0;
        if (c_map[x + y * geom.s_x] == ' ')
          continue;
        if (x < geom.s_x - 1 && c_map[x + 1 + y * geom.s_x] == ' ')
          check |= 8;
        if (x > 0 && c_map[x - 1 + y * geom.s_x] == ' ')
          check |= 2;
        if (y < geom.s_y - 1 && c_map[x + (y + 1) * geom.s_x] == ' ')
          check |= 1;
        if (y > 0 && c_map[x + (y - 1) * geom.s_x] == ' ')
          check |= 4;
        switch (check) {
        case 5:
        case 7:
        case 10:
        case 11:
        case 13:
        case 14:
        case 15:
          c_map[x + y * geom.s_x] = ' ';
          if (x > 0)
            x--;
          if (y > 0)
            y--;
          break;
        default:
          break;
        };
      }
    }
  }
  if (cnt2) {
    FILE *handle;
    char name[100];
    sprintf(name, "%s.geom", filename);
    handle = fopen(name, "w");
    fprintf(handle, "size = %i %i\n", geom.s_x, geom.s_y);
    fprintf(handle, "length = %lf %lf\n", geom.l_x, geom.l_y);
    fprintf(handle, "velocity = %lf %lf\n", geom.v_x, geom.v_y);
    fprintf(handle, "pressure = %lf\n", geom.p);
    fprintf(handle, "geometry = free\n");
    if (geom.pre == 5) {
      for (int j = geom.s_y - 1; j >= 0; --j) {
        for (int i = 0; i < geom.s_x; ++i)
          fprintf(handle, "%c", c_map[i + j * geom.s_x]);
        fprintf(handle, "\n");
      }
    } else {
      int map_cnt = 0;
      if (geom.pre == 4 && pre_data1 > 1) {
        c_map = new char[pre_data1*pre_data1];
        for (int i = 0; i < pre_data1*pre_data1; ++i)
          c_map[i] = ' ';
        int lx = -1;
        int ly = -1;
        for (int k = -pre_data1/2; k <= pre_data1/2; ++k) {
          int px = (float(pre_data1))/2.0 + float(k) * cos(pre_a) - 0.5;
          int py = (float(pre_data1))/2.0 + float(k) * sin(pre_a) - 0.5;
          if (px >= pre_data1) continue;
          if (px < 0) continue;
          if (py >= pre_data1) continue;
          if (py < 0) continue;
          c_map[px + py*pre_data1] = '#';
          if (px < pre_data1 - 1 && ly != py)
            c_map[px + py*pre_data1 + 1] = '#';
          if (py < pre_data1 - 1 && lx != px)
            c_map[px + py*pre_data1 + pre_data1] = '#';
          if (k == pre_data1/2) {
            int map = 0;
            int cnt = 0;
            if (py < pre_data1 - 1  && c_map[px + py*pre_data1 + pre_data1] != ' ') {
              map |= 1;
            }
            if (px > 0 && py < pre_data1 - 1  && c_map[px + py*pre_data1 + pre_data1 - 1] != ' ') {
              map |= 2;
            }
            if (px > 0 && c_map[px + py*pre_data1 - 1] != ' ') {
              map |= 4;
            }
            if (px > 0 && py > 0 && c_map[px + py*pre_data1 - pre_data1 - 1] != ' ') {
              map |= 8;
            }
            if (py > 0 && c_map[px + py*pre_data1 - pre_data1] != ' ') {
              map |= 16;
            }
            if (px < pre_data1 - 1 && py > 0  && c_map[px + py*pre_data1 - pre_data1 + 1] != ' ') {
              map |= 32;
            }
            if (px < pre_data1 - 1 && c_map[px + py*pre_data1 + 1] != ' ') {
              map |= 64;
            }
            if ((map & 1) && !(map & 2)) {
              if (map & 64)  { 
                c_map[px + py*pre_data1 + pre_data1 + 1] = '#'; map |= 128;
              }
              else if (map & 4)  { 
                c_map[px + py*pre_data1 + pre_data1 - 1] = '#';
                map |= 2;
              }
              else c_map[px + py*pre_data1 + pre_data1] = ' ';
            }
            if ((map & 64) && !((map & 32) || (map & 128))) {
              if (map & 16) c_map[px + py*pre_data1 - pre_data1 + 1] = '#';
              else c_map[px + py*pre_data1 + 1] = ' ';
            }
          }
          lx = px;
          ly = py;
        }
      }
      fprintf(handle, "#");
      for (int i = 0; i < geom.s_x - 2; ++i)
        switch (geom.pre) {
        case 1:
        case 2:
        case 3:
        case 4:
          fprintf(handle, "#");
          break;
        default:
          fprintf(handle, "I");
          break;
        };
      fprintf(handle, "#\n");
      for (int j = 0; j < geom.s_y - 2; ++j) {
        switch (geom.pre) {
        case 1:
          fprintf(handle, "V");
          for (int i = 0; i < geom.s_x - 2; ++i)
            fprintf(handle, " ");
          fprintf(handle, "O\n");
          break;
        case 2:
          fprintf(handle, "-");
          for (int i = 0; i < geom.s_x - 2; ++i)
            fprintf(handle, " ");
          fprintf(handle, "O\n");
          break;
        case 3:
          if (j + 1 < geom.s_y - pre_data1)
            fprintf(handle, "-");
          else
            fprintf(handle, "#");
          for (int i = 0; i < geom.s_x - 2; ++i) {
            if (j + 1 < geom.s_y - pre_data1 || i > pre_data1)
              fprintf(handle, " ");
            else
              fprintf(handle, "#");
          }
          fprintf(handle, "O\n");
          break;
        case 4:
          fprintf(handle, "-");
          for (int i = 0; i < geom.s_x - 2; ++i) {
            int up = (geom.s_y - pre_data1)/2;
            if (j + 1 < up || j + 1 >= up + pre_data1)
              fprintf(handle, " ");
            else if (i + 1 < up || i + 1 >=  up + pre_data1)
              fprintf(handle, " ");
            else
              fprintf(handle, "%c", c_map[map_cnt++]);
            /*if (j + 1 < geom.s_y / 4 || j + 1 > geom.s_y * 3 / 4)
              fprintf(handle, " ");
            else if (i == j + geom.s_y / 4 || i == j + geom.s_y / 4 + 1)
              fprintf(handle, "#");
            else if (j == geom.s_y / 4 && i + 1 == j + geom.s_y / 4)
              fprintf(handle, "#");
            else if (j + 2 == geom.s_y * 3 / 4 && i - 2 == j + geom.s_y / 4)
              fprintf(handle, "#");
            else
              fprintf(handle, " ");
            */
          }
          fprintf(handle, "O\n");
          break;
        default:
          fprintf(handle, "#");
          for (int i = 0; i < geom.s_x - 2; ++i)
            fprintf(handle, " ");
          fprintf(handle, "#\n");
          break;
        };
      }
      for (int i = 0; i < geom.s_x; ++i)
        fprintf(handle, "#");
    }
    fclose(handle);
  }

  if (c_map)
    delete[] c_map;
  if (!cnt && !cnt2) {
    printf("Greetings to you. This is a recorded announcement as I’m afraid "
           "we’re all\n out at the moment. The commercial council of Magrathea "
           "thanks you for your\n esteemed visit, but regrets that the entire "
           "planet is closed for business.\n\nThank you.\n\nIf you would like "
           "to leave your name, and a planet where you can be contacted,\n "
           "kindly speak when you hear the tone.\n");
    printf("\nWe would like to assure you that as soon as our business is "
           "resumed,\n announcements will be made in all fashionable magazines "
           "and colour supplements\n when our clients will once again be able "
           "to select from all that’s best in\n contemporary geography. "
           "Meanwhile we thank our clients for their kind interest\n and would "
           "ask them to leave. Now.\n");
    printf("\nIt is most gratifying that your enthusiasm for our planet "
           "continues unabated.\n And so we would like to assure you that the "
           "guided missiles currently\n converging with your ship are part of "
           "a special service we extend to all of\n our most enthusiastic "
           "client... And the fully armed nuclear warheads are,\n of course, "
           "merely a courtesy detail. We look forward to your custom in "
           "future\n lives. Thank you.\n");
  }

  return bye();
}
