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
    std::cout << "Reading data" << std::endl;
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
            if(count){
                std::cout << count << std::endl;
            }
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
    return;
}

void fit_peak_ge(TH1D* input_hist, double search_min, double search_max, double* mean, double* error){

    // Initialise the fit
    TF1* ge_fit = new TF1("gausexpo", "gaus[0]", search_min - 20, search_max + 20);

    // Make a first guess
    input_hist->SetAxisRange(search_min, search_max);
    double guess_mean = input_hist->GetMean();
    double guess_sigma = input_hist->GetRMS();

    input_hist->Fit("gaus", "PI", "", guess_mean - 2.0 * guess_sigma, guess_mean + 2.0 * guess_sigma);
    TF1* fit_func = input_hist->GetFunction("gaus");
    double norm = fit_func->GetParameter(0);
    guess_mean = fit_func->GetParameter(1);
    guess_sigma = fit_func->GetParameter(2);
    double guess_error = fit_func->GetParError(1);
    *mean = guess_mean;
    *error = guess_sigma;
}

void plot_channel_hist(std::string inputFile){
     // Open the file input name and check that it has opened properly
    std::ifstream input_data;
    input_data.open(inputFile.c_str());
    if(!input_data.is_open()){
        std::cout << "!!!!! Cannot find " << inputFile.c_str() << "!!!!!" << std::endl;
        return;
    }
    std::cout << "Reading data" << std::endl;
    // Create a buffer string to load data into
    std::string buffer;

    // Counter for the bin that is currently being read from the data file
    int bin = 0;
    TH1D* hist = new TH1D("", "", 4096, 0, 4096);
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
            if(count){
                std::cout << count << std::endl;
            }
            bin++;
            if (bin <= nbins) {
                hist->SetBinContent(bin, count);
            }
        }
    }

    hist->Draw();
}
