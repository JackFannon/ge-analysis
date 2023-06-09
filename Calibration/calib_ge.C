#include <string>
#include "../all_headers.h"
#include "functions.C"

// Config variables
// data_type
// 0: Old Ge data
// 1: New Ge data, no attenuator
// 3: New Ge data with attenuator (-4 db)
int data_type = 0;
bool Ni_sum_flag = false;

const int nhists = 1;

std::string colour[nhists] = {
  "Co60_1"
};

int source_type[1] = {
  1
};

std::string file_name[nhists] = {
  "../data/LINAC2021/OldGeOldDAQLIN06MeV+Co60_20210805_2306.TXT"////2021linac
};

// Rough corresponding run number in order to plot with other Ge data with LINAC beam
int runno_dummy[nhists] = {
  81586
};

TH1D * h[nhists];
TH1D * h_calib[nhists];

void calib_ge(){

  int fontid = 132;
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatFont(fontid);
  gStyle->SetLabelFont(fontid, "XYZ");
  gStyle->SetLabelFont(fontid, "");
  gStyle->SetTitleFont(fontid, "XYZ");
  gStyle->SetTitleFont(fontid, "");
  gStyle->SetTitleOffset(1.2, "XYZ");
  gStyle->SetTextFont(fontid);


}
