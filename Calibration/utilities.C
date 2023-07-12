#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
//#include "../all_headers.h"

void read_data_into_hist(std::string inputname, TH1D* hist){
     // Open the file input name and check that it has opened properly
    std::ifstream input_data;
    input_data.open(inputname.c_str());
    if(!input_data.is_open()){
        std::cout << "!!!!! Cannot find " << inputname.c_str() << "!!!!!" << std::endl;
        return;
    }
    // Create a buffer string to load data into
    std::string buffer;

    // Counter for the bin that is currently being read from the data file
    int bin = 0;
    // Max number of bins that we should read to
    int nbins = hist->GetNbinsX();

    // Need to skip the first 9 lines of the data file as these are all meta-data comments
    for (int i = 0; i < 9; i++) {
        std::getline(input_data, buffer);
    }

    // String to store the channel number in (this is not used, but is good to be named something other than "useless variable")
    int channel;

    // Format of data is as follows:
    // Channel, RSV1, RSV2, ..., RSV9, RSV10
    // We're not interested in "Channel" but are interested in RSV1-10

    // Loop over the rest of the file
    while (!input_data.eof()) {
        std::getline(input_data, buffer);

        // Put contents on the line into a stringstream
        std::stringstream line(buffer);

        std::string token;
        // Load the channel number from "line" into "channel"

        line >> channel;
        // Loop over the remaining RSV numbers in the line
        for (int i = 0; i < 10; i++) {
            int count = -9999;
            if(line.peek() == ','){
                line.ignore();
            }
            line >> count;
            bin++;
            if (bin <= nbins) {
                hist->SetBinContent(bin, count);
            }
        }
    }
}

void set_style(int fontid){
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
    gStyle->SetOptFit();
    gStyle->GetPadTopMargin();
    return;
}

void fit_peak_ge(TH1D* input_hist, double search_min, double search_max, double* mean, double* error){

    // Initialise the fit
    TF1* ge_fit = new TF1("gauslin", "gaus(0) + pol1(3)", search_min, search_max);
    ge_fit->SetParLimits(0, 0., pow(10., 6));
    ge_fit->SetParLimits(2, .1, 10.);
    ge_fit->SetParameter(0, 100.);

    // Make a first guess
    input_hist->SetAxisRange(search_min, search_max);
    double guess_mean = input_hist->GetMean();
    double guess_sigma = input_hist->GetRMS();

    ge_fit->SetParLimits(1, guess_mean - 3, guess_mean + 3);
    ge_fit->SetParameter(1, guess_mean);
    ge_fit->SetParameter(2, guess_sigma);

    for (int i = 0; i < 10; i++) {
        input_hist->Fit( "gauslin", "LIRQ");
    }

    input_hist->Fit( "gauslin", "LIR");
    double norm = ge_fit->GetParameter(0);
    guess_mean = ge_fit->GetParameter(1);
    guess_sigma = ge_fit->GetParameter(2);
    double guess_error = ge_fit->GetParError(1);
    *mean = guess_mean;
    *error = guess_sigma;
}

TH1D* plot_channel_hist(std::string inputFile, std::string directory){
    // Open the file input name and check that it has opened properly
    std::ifstream input_data;
    input_data.open(directory + inputFile.c_str());
    if(!input_data.is_open()){
        std::cout << "!!!!! Cannot find " << inputFile.c_str() << "!!!!!" << std::endl;
    }
    // Create a buffer string to load data into
    std::string buffer;

    // Counter for the bin that is currently being read from the data file
    int bin = 0;
    TH1D* hist = new TH1D(inputFile.c_str(), inputFile.c_str(), 4096, 0, 4096);
    // Max number of bins that we should read to
    int nbins = hist->GetNbinsX();

    // Need to skip the first 9 lines of the data file as these are all meta-data comments
    for (int i = 0; i < 9; i++) {
        std::getline(input_data, buffer);
    }

    // String to store the channel number in (this is not used, but is good to be named something other than "useless variable")
    int channel;

    // Format of data is as follows:
    // Channel, RSV1, RSV2, ..., RSV9, RSV10
    // We're not interested in "Channel" but are interested in RSV1-10

    // Loop over the rest of the file
    while (!input_data.eof()) {
        std::getline(input_data, buffer);

        // Put contents on the line into a stringstream
        std::stringstream line(buffer);

        std::string token;
        // Load the channel number from "line" into "channel"

        line >> channel;
        // Loop over the remaining RSV numbers in the line
        for (int i = 0; i < 10; i++) {
            int count = -9999;
            if(line.peek() == ','){
                line.ignore();
            }
            line >> count;
            bin++;
            if (bin <= nbins) {
                hist->SetBinContent(bin, count);
            }
        }
    }

    return hist;
}


std::vector<std::string> load_data(std::string filenames, std::string directory){

    std::vector<std::string> file_list;

    std::ifstream input_file;
    input_file.open(directory + filenames);

    if(!input_file.is_open()){
        std::cerr << "Could not open " << directory + filenames << std::endl;
    }

    std::string buffer;
    while(!input_file.eof()){
        std::getline(input_file, buffer);
        std::stringstream input_line(buffer);
        std::string file_name;
        input_line >> file_name;
        if(!file_name.empty()){
            file_list.push_back(file_name);
        }
    }

    return file_list;
}
