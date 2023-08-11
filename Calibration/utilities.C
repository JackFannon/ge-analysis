#include "TH1.h"
#include "TStyle.h"
#include "TF1.h"
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

// Reads a raw data file, inputname, and populates the TH1F, hist, with
//   channel counts. This is not calibrated, just channel vs counts.
void read_data_into_hist(std::string inputname, TH1F* hist){
     // Open the file input name and check that it has opened properly
    std::ifstream input_data;
    input_data.open(inputname.c_str());
    if(!input_data.is_open()){
        std::cout << "!!!!! Cannot find " << inputname.c_str() << "!!!!!" << std::endl;
        return;
    }
    // Create a buffer string to load data into
    std::string buffer;

    // Need to check what format the data is in - new detector vs old detector
    std::getline(input_data, buffer);

    std::istringstream line(buffer);

    std::string firstWord;

    std::getline(line, firstWord, ',');

    std::cout << firstWord << std::endl;


    try{
        if(firstWord != "1SPECTRUM" && firstWord != "1Channel"){
            throw firstWord;
        }
    }
    catch (std::string &word){
        std::cout << "Error: Data type may not be supported. Check the documentation for supported types, functionality may need to be added to read your data." << word << std::endl;
    }
    catch (const char* word){
        std::cout << "Error: Data type may not be supported. Check the documentation for supported types, functionality may need to be added to read your data." << word << std::endl;
    }




    if(firstWord == "SPECTRUM"){
        // WE HAVE DATA FROM THE NEW DETECTOR
        std::cout << "NEW DETECTOR" << std::endl;
    }
    else if(firstWord == "1Channel"){
        // WE HAVE DATA FROM THE OLD DETECTOR
        std::cout << "OLD DETECTOR" << std::endl;
    }

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
            if (bin <= nbins){
                hist->SetBinContent(bin, count);
            }
        }
    }
}

// Sets the drawing style
// NOTE - Not sure if any of this actually works. ROOT seems to ignore gStyle->SetOptStat(0) anyway
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

// Attempts to fit a Gaussian + linear fit between search_min and search_max
// Populates the variables mean and error with the mean value and standard deviation of the
//   Gaussian fit.
void fit_peak_ge(TH1F* input_hist, double search_min, double search_max, double* mean, double* error){

    // Initialise the fit
    TF1* ge_fit = new TF1("gauslin", "gaus(0) + pol1(3)", search_min, search_max);
    ge_fit->SetParLimits(0, 0., pow(10., 6));
    ge_fit->SetParLimits(2, .1, 10.);
    ge_fit->SetParameter(0, 100.);
    //ge_fit->SetRange(search_min, search_max);
    std::cout << __LINE__ << std::endl;

    TH1F* copy_hist = (TH1F*)input_hist->Clone("copy");
    copy_hist->SetAxisRange(0.98 * search_min, 1.02 * search_max);
    // Make a first guess
    double guess_mean = copy_hist->GetMean();
    double guess_sigma = copy_hist->GetRMS();

    std::cout << guess_mean << std::endl;

    ge_fit->SetParLimits(1, guess_mean, guess_mean);
    ge_fit->SetParameter(1, guess_mean);
    ge_fit->SetParameter(2, guess_sigma);

    for (int i = 0; i < 10; i++){
        copy_hist->Fit("gauslin", "LIRQ");
    }

    copy_hist->Fit("gauslin", "LIR");
    double norm = ge_fit->GetParameter(0);
    guess_mean = ge_fit->GetParameter(1);
    guess_sigma = ge_fit->GetParameter(2);
    double guess_error = ge_fit->GetParError(1);
    *mean = guess_mean;
    *error = guess_sigma;
    delete copy_hist;
}

// Plots a histogram of channel number vs counts
TH1F* plot_channel_hist(std::string inputFile, std::string directory){
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
    TH1F* hist = new TH1F(inputFile.c_str(), inputFile.c_str(), 4096, 0, 4096);
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

// Returns a vector containing a list of paths for the raw Ge data
// filename - Filename of the list of files to read
// directory - location of the file with name filename
std::vector<std::string> load_data(std::string filename, std::string directory){

    std::vector<std::string> file_list;

    std::ifstream input_file;
    input_file.open((directory + filename).c_str());

    if(!input_file.is_open()){
        std::cerr << "Could not open " << directory + filename << std::endl;
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
