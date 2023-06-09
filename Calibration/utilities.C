#include <sstream>
#include <string>
#include <fstream>
#include "../all_headers.h"

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
    std::string channel;

    // Format of data is as follows:
    // Channel, RSV1, RSV2, ..., RSV9, RSV10
    // We're not interested in "Channel" but are interested in RSV1-10

    // Loop over the rest of the file
    while (!input_data.eof()) {
        std::getline(input_data, buffer);

        // Put contents on the line into a stringstream
        std::stringstream line(buffer);

        // Load the channel number from "line" into "channel"
        line >> channel;

        // Loop over the remaining RSV numbers in the line
        for (int i = 0; i < 10; i++) {
            int count = 0;
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
    return;
}
