
#include "string"
#include "iostream"

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TH2.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TFitResult.h"
#include "TMinuit.h"

using namespace std;


class Detection {
    public:
        UInt_t detector_id;
        Int_t chip_id;
        Float_t momentum;
        Float_t pt;
        Float_t path_length;
        Float_t deltaE;
        Float_t eta;
        Float_t kaon;
        Float_t clean_RP_conf;
        Detection() {  
        }
};

class Track {
    public:
        vector<Detection*> detections;
        Float_t p_t;
        Float_t total_deltaE = 0;
        Float_t total_path_length = 0;
        Float_t total_p = 0;
        Float_t average_p;
        Float_t dEdx;

        Track() {  
        }

        void reset() {
            detections.clear();
        }

        void calculate_all() {
            calculate_total_deltaE();
            calculate_total_path_length();
            calculate_total_p();
            calculate_average_p();
            calculate_dEdx();
        }

        void calculate_total_deltaE() {
            Float_t total = 0;
            for (Detection* detection : detections) {
                total += detection->deltaE;
            }
            total_deltaE = total;
        }

        void calculate_total_path_length() {
            Float_t total = 0;
            for (Detection* detection : detections) {
                total += detection->path_length;
            }
            total_path_length = total;
        }

        void calculate_total_p() {
            Float_t total = 0;
            for (Detection* detection : detections) {
                total += detection->momentum;
            }
            total_p = total;
        }

        void calculate_average_p() {
            if (total_p == 0) {
                calculate_total_p();
            }
            average_p = total_p/detections.size();
        }

        void calculate_dEdx(){
            float sum = 0;
            for (Detection* detection : detections) {
                sum += (detection->deltaE)/(detection->path_length);
            }
            dEdx = sum/detections.size();
        }
};

string set_length(float value, int length) {
    string value_string = to_string(value);
    int value_string_length = value_string.length();
    int difference = length - value_string_length;
    string final_string = "";
    if (difference > 0) {
        for (int i=0; i<difference; ++i) {
            final_string += " ";
        }
        final_string += value_string;
    } else if (difference < 0) {
        for (int i=0; i<value_string_length+difference; i++) {
            final_string += value_string[i];
        }
    } else {
        final_string += value_string;
    }
    return final_string;
}

string set_length(int value, int length) {
    string value_string = to_string(value);
    int difference = length - value_string.length();
    if (difference < 0) {
        throw invalid_argument("length is too short or value is too long"); 
    }
    string final_string = "";
    for (int i=0; i<difference; ++i) {
        final_string += "0";
    }
    final_string += value_string;
    return final_string;
}

int main() {

    const string filenames[] = {"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_00.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_01.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_02.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_03.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_04.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_05.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_06.root"};
    
/*
{"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_00.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_01.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_02.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_03.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_04.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_05.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_06.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_07.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_08.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_09.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_10.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_11.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_12.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_13.root", "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_14.root"}
*/

    TH2F* dEdx = new TH2F("dEdx", ";p;dEdx",200,0,2.5,200,0,120);
    auto path_vs_dE = new TH2F("path_vs_dE", ";path length;delta E",200,0,0.01,200,0,0.5);
    auto track_dEdx = new TH1F("track_dEdx", ";dE/dx;detections",200,0,100);
    
    //map<int, int> detectors;


    for (string filename : filenames) {
        
        cout << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("mytree");
        TTreeReader Reader(tree);

        //gStyle->SetOptStat(0);
        //gStyle->SetPalette(kCividis);


        TTreeReaderValue<UInt_t> tree_detector_id(Reader, "tree_detector_id");
        TTreeReaderValue<Int_t> tree_chip_id(Reader, "tree_chip_id");
        TTreeReaderValue<Float_t> tree_momentum(Reader, "tree_momentum");
        TTreeReaderValue<Float_t> tree_pt(Reader, "tree_pt");
        TTreeReaderValue<Float_t> tree_path_length(Reader, "tree_path_length");
        TTreeReaderValue<Float_t> tree_deltaE(Reader, "tree_deltaE");
        TTreeReaderValue<Float_t> tree_eta(Reader, "tree_eta");
        TTreeReaderValue<Bool_t> tree_kaon(Reader, "tree_kaon");
        TTreeReaderValue<Bool_t> tree_clean_RP_conf(Reader, "tree_clean_RP_conf");

        float previous_tree_pt = 0;

        Track* track = new Track();


        while (Reader.Next()) {
            if (*tree_pt != previous_tree_pt && previous_tree_pt != 0) {

                track->calculate_all();
                dEdx->Fill(track->average_p, track->dEdx);
                track->reset();
                
            }
            previous_tree_pt = *tree_pt;
            

            Detection* detection = new Detection();
            detection->detector_id = *tree_detector_id;
            detection->chip_id = *tree_chip_id;
            detection->momentum = *tree_momentum;
            detection->pt = *tree_pt;
            detection->path_length = 0.1*(*tree_path_length);
            detection->deltaE = *tree_deltaE;
            detection->eta = *tree_eta;
            detection->kaon = *tree_kaon;
            detection->clean_RP_conf = *tree_clean_RP_conf;

            track->detections.push_back(detection);

        }
    }



    auto base_canvas = new TCanvas("base_canvas","base_canvas");

    dEdx->Draw("Colz");
    
    base_canvas->Print("dEdx_graph.pdf");

}
